// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"

#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChElementCableANCF.h"
#include "chrono/fea/ChBuilderBeam.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include <chrono_interiorpoint/ChInteriorPoint.h>

using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

//#define BUILD_CABLES
#define BUILD_TETRAHEDRONS

int main(int argc, char* argv[]) {

    // Create a Chrono::Engine physical system
    ChSystemNSC my_system;


    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"IP FEA contacts", core::dimension2d<u32>(800, 600), false, true, true, irr::video::EDT_OPENGL);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(0, (f32)0.6, -1));
    application.AddLightWithShadow(core::vector3df(1.5, 5.5, -2.5), core::vector3df(0, 0, 0), 3, 2.2, 7.2, 40, 512,
        video::SColorf(1, 1, 1));

	application.SetContactsDrawMode(ChIrrTools::CONTACT_NORMALS);
	application.SetSymbolscale(0.25);

    //
    // CREATE THE PHYSICAL SYSTEM
    //


    //collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0); // not needed, already 0 when using ChSystemSMC
	collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
	collision::ChCollisionModel::SetDefaultSuggestedMargin(0.001);

    // Use this value for an outward additional layer around meshes, that can improve
    // robustness of mesh-mesh collision detection (at the cost of having unnatural inflate effect)
    double sphere_swept_thickness = 0.002;

    // Create the surface material, containing information
    // about friction etc. 
    // It is a NSC (penalty) material that we will assign to 
    // all surfaces that might generate contacts.

    auto mysurfmaterial = std::make_shared<ChMaterialSurfaceNSC>();
    mysurfmaterial->SetFriction(0.4f);
    mysurfmaterial->SetCompliance(0.0000005f);
    mysurfmaterial->SetComplianceT(0.0000005f);
    mysurfmaterial->SetDampingF(0.2f);

    // Create a floor as a simple collision primitive:
    auto mfloor = std::make_shared<ChBodyEasyBox>(2, 0.1, 2, 2700, true);
    mfloor->SetPos(ChVector<>(0,-0.1/2, 0));
    mfloor->SetBodyFixed(true);
    mfloor->SetMaterialSurface(mysurfmaterial);
    my_system.Add(mfloor);

    auto masset_texture = std::make_shared<ChTexture>();
    masset_texture->SetTextureFilename(GetChronoDataFile("concrete.jpg"));
    mfloor->AddAsset(masset_texture);



    // two falling objects:

    //auto mcube = std::make_shared<ChBodyEasyBox>(0.1, 0.1, 0.1, 2700, true);
    //mcube->SetPos(ChVector<>(0.6, 0.5, 0.6));
    //mcube->SetMaterialSurface(mysurfmaterial);
    //my_system.Add(mcube);

    //auto msphere = std::make_shared<ChBodyEasySphere>(0.1, 2700, true);
    //msphere->SetPos(ChVector<>(0.8, 0.5, 0.6));
    //msphere->SetMaterialSurface(mysurfmaterial);
    //my_system.Add(msphere);

#ifdef BUILD_TETRAHEDRONS
    //
    // Example 1: tetrahedrons, with collisions
    // 

    // Create a mesh. We will use it for tetahedrons.

    auto my_mesh = std::make_shared<ChMesh>();

    // 1) a FEA tetahedron(s):

    // Create a material, that must be assigned to each solid element in the mesh,
    // and set its parameters
    auto mmaterial = std::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(0.01e9);  // rubber 0.01e9, steel 200e9
    mmaterial->Set_v(0.3);
    mmaterial->Set_RayleighDampingK(0.003);
    mmaterial->Set_density(1000);

    if (false) {
        for (int k = 0; k<3; ++k)
            for (int j = 0; j<3; ++j)
                for (int i = 0; i<3; ++i) {
                    // Creates the nodes for the tetahedron
                    ChVector<> offset(j*0.21, i*0.21, k*0.21);
                    auto mnode1 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0.1, 0) + offset);
                    auto mnode2 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0.1, 0.2) + offset);
                    auto mnode3 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0, 0.3, 0) + offset);
                    auto mnode4 = std::make_shared<ChNodeFEAxyz>(ChVector<>(0.2, 0.1, 0) + offset);

                    my_mesh->AddNode(mnode1);
                    my_mesh->AddNode(mnode2);
                    my_mesh->AddNode(mnode3);
                    my_mesh->AddNode(mnode4);

                    auto melement1 = std::make_shared<ChElementTetra_4>();
                    melement1->SetNodes(mnode1,
                        mnode2,
                        mnode3,
                        mnode4);
                    melement1->SetMaterial(mmaterial);

                    my_mesh->AddElement(melement1);
                }
    }

    if (true) {
        for (int i = 0; i<4; ++i) {
            try
            {
                ChCoordsys<> cdown(ChVector<>(0, -0.4, 0));
                ChCoordsys<> crot(VNULL, Q_from_AngAxis(CH_C_2PI * ChRandom(), VECT_Y) * Q_from_AngAxis(CH_C_PI_2, VECT_X));
                ChCoordsys<> cydisp(ChVector<>(-0.3, 0.1 + i*0.1, -0.3));
                ChCoordsys<> ctot = cdown >> crot >> cydisp;
                ChMatrix33<> mrot(ctot.rot);
                ChMeshFileLoader::FromTetGenFile(my_mesh, GetChronoDataFile("fea/beam.node").c_str(), GetChronoDataFile("fea/beam.ele").c_str(), mmaterial, ctot.pos, mrot);
            }
            catch (ChException myerr) {
                GetLog() << myerr.what();
                return 0;
            }
        }
    }


    //Create the contact surface(s). 
    //In this case it is a ChContactSurfaceMesh, that allows mesh-mesh collsions.

    auto mcontactsurf = std::make_shared<ChContactSurfaceMesh>();
    my_mesh->AddContactSurface(mcontactsurf);
    mcontactsurf->AddFacesFromBoundary(sphere_swept_thickness); // do this after my_mesh->AddContactSurface
    mcontactsurf->SetMaterialSurface(mysurfmaterial); // use the NSC penalty contacts
    
    // Remember to add the mesh to the system!

    //auto mcontactclouda = std::make_shared<ChContactSurfaceNodeCloud>();
    //my_mesh->AddContactSurface(mcontactclouda);
    //mcontactclouda->AddAllNodes(0.005); // use larger point size to match beam section radius
    //mcontactclouda->SetMaterialSurface(mysurfmaterial);
    //

    my_system.Add(my_mesh);

	// ==Asset== attach a visualization of the FEM mesh.
	// This will automatically update a triangle mesh (a ChTriangleMeshShape
	// asset that is internally managed) by setting  proper
	// coordinates and vertex colours as in the FEM elements.
	// Such triangle mesh can be rendered by Irrlicht or POVray or whatever
	// postprocessor that can handle a coloured ChTriangleMeshShape).
	// Do not forget AddAsset() at the end!

	auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
	mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
	mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
	mvisualizemesh->SetSmoothFaces(true);
	my_mesh->AddAsset(mvisualizemesh);

	auto mvisualizemeshcoll = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
	mvisualizemeshcoll->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_CONTACTSURFACES);
	mvisualizemeshcoll->SetWireframe(true);
	mvisualizemeshcoll->SetDefaultMeshColor(ChColor(1, 0.5, 0));
	my_mesh->AddAsset(mvisualizemeshcoll);


#endif

#ifdef BUILD_CABLES
    //
    // Example 2: beams, with collisions
    // 

    // Create a mesh. We will use it for beams only. 

    auto my_mesh_beams = std::make_shared<ChMesh>();

    // 2) an ANCF cable:

    auto msection_cable2 = std::make_shared<ChBeamSectionCable>();
    msection_cable2->SetDiameter(0.05);
    msection_cable2->SetYoungModulus(0.01e9);
    msection_cable2->SetBeamRaleyghDamping(0.05);

    ChBuilderBeamANCF builder;

    double cable_height = 0.0;
    for (auto cable_sel = 0; cable_sel < 10; ++cable_sel)
    {
        cable_height += 0.1;

        builder.BuildBeam(my_mesh_beams,	// the mesh where to put the created nodes and elements 
            msection_cable2,// the ChBeamSectionCable to use for the ChElementCableANCF elements
            10,				// the number of ChElementCableANCF to create
            ChVector<>(-0.25, cable_height, -0.25),		// the 'A' point in space (beginning of beam)
            ChVector<>(0.25, cable_height, 0.25));	// the 'B' point in space (end of beam)

        cable_height += 0.1;

        builder.BuildBeam(my_mesh_beams,	// the mesh where to put the created nodes and elements 
            msection_cable2,// the ChBeamSectionCable to use for the ChElementCableANCF elements
            10,				// the number of ChElementCableANCF to create
            ChVector<>(-0.25, cable_height, 0.25),		// the 'A' point in space (beginning of beam)
            ChVector<>(0.25, cable_height, -0.25));	// the 'B' point in space (end of beam)

        std::cout << cable_height << std::endl;
    }

	// ==Asset== attach a visualization of the FEM mesh.
	// This will automatically update a triangle mesh (a ChTriangleMeshShape
	// asset that is internally managed) by setting  proper
	// coordinates and vertex colours as in the FEM elements.
	// Such triangle mesh can be rendered by Irrlicht or POVray or whatever
	// postprocessor that can handle a coloured ChTriangleMeshShape).
	// Do not forget AddAsset() at the end!

    // Create the contact surface(s). 
    // In this case it is a ChContactSurfaceNodeCloud, so just pass 
    // all nodes to it.
    auto mcontactcloud = std::make_shared<ChContactSurfaceNodeCloud>();
    my_mesh_beams->AddContactSurface(mcontactcloud);
    mcontactcloud->AddAllNodes(0.001); // use larger point size to match beam section radius
    mcontactcloud->SetMaterialSurface(mysurfmaterial);

    auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
    mvisualizemesh->SetSmoothFaces(true);
    my_mesh_beams->AddAsset(mvisualizemesh);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh_beams);

	auto mvisualizemeshbeam = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
	mvisualizemeshbeam->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
	mvisualizemeshbeam->SetColorscaleMinMax(0.0, 5.50);
	mvisualizemeshbeam->SetSmoothFaces(true);
	my_mesh->AddAsset(mvisualizemeshbeam);

	auto mvisualizemeshbeamnodes = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh_beams.get()));
	mvisualizemeshbeamnodes->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
	mvisualizemeshbeamnodes->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
	mvisualizemeshbeamnodes->SetSymbolsThickness(0.008);
	my_mesh->AddAsset(mvisualizemeshbeamnodes);

#endif

    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
    // in the system. These ChIrrNodeAsset assets are 'proxies' to the Irrlicht meshes.
    // If you need a finer control on which item really needs a visualization proxy in
    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.
    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!

    application.AssetUpdateAll();

    // Use shadows in realtime view
    application.AddShadowAll();


    // Mark completion of system construction
    my_system.SetupInitial();


    //
    // THE SOFT-REAL-TIME CYCLE
    //
    my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

    // Change solver to IP
    auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
    auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
    my_system.SetStabSolver(ip_solver_stab);
    my_system.SetSolver(ip_solver_speed);
    ip_solver_speed->RecordHistory(false);
    ip_solver_speed->SetNullPivotDetection(true, 1e-18);
    ip_solver_speed->SetVerbose(true);

    application.GetSystem()->Update();
    application.SetPaused(true);

    application.SetTimestep(0.01);

    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        application.DoStep();

        application.EndScene();
    }

    return 0;
}
