//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be 
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

//
//   Demo code about 
//
//     - FEA for 3D beams



// Include some headers used by this tutorial...

#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/solver/ChSolverMINRES.h"

#include "chrono_fea/ChElementBeamEuler.h"
#include "chrono_fea/ChBuilderBeam.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include "core/ChCSR3Matrix.h"
#include "chrono_modelreduction/ChEigenAnalysis.h"
#include <typeinfo>

//#include "chrono_matlab/ChMatlabEngine.h"
//#include "chrono_matlab/ChSolverMatlab.h"

// Remember to use the namespace 'chrono' because all classes 
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;



int main(int argc, char* argv[])
{
    // Create a Chrono::Engine physical system
    ChSystem my_system;


    // Create the Irrlicht visualization (open the Irrlicht device, 
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Beam modes", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(-0.1f, 0.2f, -0.2f));



    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();


    // Create a section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.

    auto msection = std::make_shared<ChBeamSectionAdvanced>();

    double beam_wy = 0.012;
    double beam_wz = 0.025;
    msection->SetAsRectangularSection(beam_wy, beam_wz);
    msection->SetYoungModulus(0.01e9);
    msection->SetGshearModulus(0.01e9 * 0.3);
    msection->SetBeamRaleyghDamping(0.000);
    msection->SetDensity(7500);
    //msection->SetCentroid(0,0.02); 
    //msection->SetShearCenter(0,0.1); 
    //msection->SetSectionRotation(45*CH_C_RAD_TO_DEG);

    //
    // Add some EULER-BERNOULLI BEAMS (the fast way!)
    //

    double beam_L = 0.2;

    ChBuilderBeam builder;

    builder.BuildBeam(my_mesh,		// the mesh where to put the created nodes and elements 
                      msection,		// the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
                      5,				// the number of ChElementBeamEuler to create
                      ChVector<>(0, 0, -0.1),		// the 'A' point in space (beginning of beam)
                      ChVector<>(beam_L, 0, -0.1),	// the 'B' point in space (end of beam)
                      ChVector<>(0, 1, 0));			// the 'Y' up direction of the section for the beam

                                                    // After having used BuildBeam(), you can retrieve the nodes used for the beam,
                                                    // For example say you want to fix the A end and apply a force to the B end:
    builder.GetLastBeamNodes().back()->SetFixed(true);
    builder.GetLastBeamNodes().front()->SetForce(ChVector<>(0, -1, 0));

    //
    // Final touches..
    // 


    // We do not want gravity effect on FEA elements in this demo
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);




    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChTriangleMeshShape
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    // postprocessor that can handle a coloured ChTriangleMeshShape).
    // Do not forget AddAsset() at the end!


    /*
    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizebeamA->SetSmoothFaces(true);
    my_mesh->AddAsset(mvisualizebeamA);
    */

    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
    mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddAsset(mvisualizebeamA);

    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddAsset(mvisualizebeamC);


    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
    // in the system. These ChIrrNodeAsset assets are 'proxies' to the Irrlicht meshes.
    // If you need a finer control on which item really needs a visualization proxy in 
    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.

    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!

    application.AssetUpdateAll();

    // Mark completion of system construction
    my_system.SetupInitial();


    // 
    // THE SOFT-REAL-TIME CYCLE
    //
    my_system.SetSolverType(ChSystem::SOLVER_MINRES);
    my_system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(460);
    my_system.SetMaxItersSolverStab(460);
    my_system.SetTolForce(1e-13);
    auto msolver = static_cast<ChSolverMINRES*>(my_system.GetSolverSpeed());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);

    application.SetTimestep(0.01);



    GetLog() << "\n\n===========EIGENPROBLEM======== \n";

    application.DoStep();

    ChEigenAnalysis eig_analysis(application);
    eig_analysis.EigenAnalysis();

    std::function<void()> simadvance_fun = std::bind(&ChEigenAnalysis::UpdateEigenMode, &eig_analysis, 1);
    //std::function<void()> render_fun = std::bind(&ChIrrApp::DrawAll, application);
    std::function<void(double)> render_fun = [&application](double dummy) { application.DrawAll(); };
    std::function<void()> pre_sim_fun = std::bind(&ChIrrApp::BeginScene, &application, true, true, irr::video::SColor(255, 0, 0, 0));

    //std::function<void(double)> render_fun = [&application, &eig_analysis](double ratio) { eig_analysis.UpdateEigenMode(ratio); application.DrawAll(); };
    ChRealtimeDualStepTimer drealtime = { application, simadvance_fun, render_fun };
    drealtime.SetPreAdvanceSimulationFunction(pre_sim_fun);


    while (application.GetDevice()->run())
    {

        //application.DrawAll();

        //eig_analysis.UpdateEigenMode();
        drealtime.DoRealtimeStep();

        application.EndScene();
    }


}


