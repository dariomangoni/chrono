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
// Authors: Alessandro Tasora
// =============================================================================
//
// FEA for 3D beams
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/collision/ChCCollisionSystemBullet.h"

#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono_irrlicht/ChIrrApp.h"

#include "chrono_interiorpoint/ChInteriorPoint.h"
#include <iomanip>

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

bool include_plasticity = false;

#define USE_NSC

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

#ifdef USE_NSC
    bool contact_model_NSC = true;
#else
    bool contact_model_NSC = false;
#endif

#ifdef USE_NSC
    ChSystemNSC my_system;
    typedef ChMaterialSurfaceNSC material_surface_class;
#else
    ChSystemSMC my_system;
    typedef ChMaterialSurfaceSMC material_surface_class;
#endif

    auto mat_surface_type = contact_model_NSC ? ChMaterialSurface::NSC : ChMaterialSurface::SMC;

    // Here set the inward-outward margins for collision shapes: should make sense in the scale of the model
    collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.001);
    collision::ChCollisionModel::SetDefaultSuggestedMargin(0.002);
    collision::ChCollisionSystemBullet::SetContactBreakingThreshold(0.0001);


    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.

    auto my_mesh = std::make_shared<ChMesh>();
    my_system.Add(my_mesh);

    // Create a section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.

    double wire_diameter = 0.010;

    auto melasticity = std::make_shared<ChElasticityCosseratSimple>();
    melasticity->SetYoungModulus(0.5e9);
    melasticity->SetGshearModulus(0.5e9 * 0.7);

    auto mdamping = std::make_shared<ChDampingCosseratLinear>();
    mdamping->SetDampingCoefficientsRe((1e-3)*ChVector<>(1, 1, 1));
    mdamping->SetDampingCoefficientsRk((1e-4)*ChVector<>(1, 1, 1));

    auto mplasticity = std::make_shared<ChPlasticityCosseratLumped>();
    mplasticity->n_yeld_Mx = std::make_shared<ChFunction_Ramp>(1, 0.01);
    mplasticity->n_yeld_My = std::make_shared<ChFunction_Ramp>(0.2, 0.001);
    mplasticity->n_yeld_Mz = std::make_shared<ChFunction_Ramp>(0.2, 0.001);

    auto msection = std::make_shared<ChBeamSectionCosserat>(melasticity, include_plasticity ? mplasticity : nullptr, mdamping);
    msection->SetDensity(1000);
    msection->SetAsCircularSection(wire_diameter);

    auto mysurfmaterial = std::make_shared<material_surface_class>();
    mysurfmaterial->SetRestitution(0.1f);
    mysurfmaterial->SetFriction(0.2f);

#ifdef USE_NSC
    mysurfmaterial->SetCompliance(0.0000005f);
    mysurfmaterial->SetComplianceT(0.0000005f);
    mysurfmaterial->SetDampingF(0.2f);
#else
    mysurfmaterial->SetYoungModulus(6e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);
#endif

    
    ///////////////
    // BEAM BUNDLE
    ///////////////
    const int beam_num = 8;
    const double bundle_radius = 0.025;
    const double bundle_length0 = 1.0;
    const double flange_thickness = 0.05;

    ChBuilderBeamIGA beam_build;
    auto fixedAsset = std::make_shared<ChColorAsset>(1.0, 0.0, 0.0); // color for fixed objects
    auto movingAsset = std::make_shared<ChColorAsset>(0.0, 0.0, 1.0); // color for moving objects

    auto mFlangeStart = std::make_shared<ChBodyEasyCylinder>((bundle_radius + wire_diameter)*1.1, flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeStart->SetPos(-0.5*VECT_Z*flange_thickness);
    mFlangeStart->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
    mFlangeStart->SetBodyFixed(true);
    mFlangeStart->AddAsset(fixedAsset);
    my_system.Add(mFlangeStart);

    auto mFlangeEnd = std::make_shared<ChBodyEasyCylinder>((bundle_radius+ wire_diameter)*1.1, flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeEnd->SetPos(VECT_Z*(bundle_length0+0.5*flange_thickness));
    mFlangeEnd->SetRot(Q_from_AngAxis(CH_C_PI_2, VECT_X));
    mFlangeEnd->AddAsset(movingAsset);
    my_system.Add(mFlangeEnd);

    auto mFlangeEndFixed = std::make_shared<ChBodyEasyBox>(2.5*(bundle_radius + wire_diameter), 2.5*(bundle_radius + wire_diameter), flange_thickness, 1000, false, true, mat_surface_type);
    mFlangeEndFixed->SetPos(mFlangeEnd->GetPos() + VECT_Z* flange_thickness);
    mFlangeEndFixed->SetBodyFixed(true);
    mFlangeEndFixed->AddAsset(fixedAsset);
    my_system.Add(mFlangeEndFixed);

    auto rotMot = std::make_shared<ChLinkMotorRotationSpeed>();
    auto speedFun = std::make_shared<ChFunction_Const>(1);
    rotMot->SetSpeedFunction(speedFun);
    //rotMot->SetAngleOffset(30.0*CH_C_DEG_TO_RAD); // not working as expected
    rotMot->Initialize(mFlangeEndFixed, mFlangeEnd, ChFrame<>(mFlangeEnd->GetPos(), QUNIT));
    my_system.Add(rotMot);
    
    // IGA beams between flanges
    for (auto beam_sel = 0; beam_sel<beam_num; ++beam_sel)
    {
        // create beams at given position
        ChVector<> beam_pos_start(bundle_radius*cos(beam_sel*CH_C_2PI / beam_num), bundle_radius*sin(beam_sel*CH_C_2PI / beam_num), 0);
        ChVector<> beam_pos_end(bundle_radius*cos(beam_sel*CH_C_2PI / beam_num), bundle_radius*sin(beam_sel*CH_C_2PI / beam_num), bundle_length0);
        beam_build.BuildBeam(my_mesh, msection, 10, beam_pos_start, beam_pos_end, VECT_Y, 3);
        auto firstNode = beam_build.GetLastBeamNodes().front();
        auto endNode = beam_build.GetLastBeamNodes().back();

        // fix beam to start flange
        auto startConstr = std::make_shared<ChLinkMateGeneric>();
        startConstr->Initialize(firstNode, mFlangeStart, false, firstNode->Frame(), firstNode->Frame());
        startConstr->SetConstrainedCoords(true, true, true, true, true, true);
        my_system.Add(startConstr);

        // fix beam to end flange
        auto endConstr = std::make_shared<ChLinkMateGeneric>();
        endConstr->Initialize(endNode, mFlangeEnd, false, endNode->Frame(), endNode->Frame());
        endConstr->SetConstrainedCoords(true, true, true, true, true, true);
        my_system.Add(endConstr);

    }

    my_mesh->SetAutomaticGravity(false);


    // VISUALIZATION
    // beams
    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
    mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddAsset(mvisualizebeamA);

    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NONE);
    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddAsset(mvisualizebeamC);


    // Irrlicht application
    ChIrrApp application(&my_system, L"Torquing IGA beams", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    auto targ = core::vector3df(mFlangeEnd->GetPos().x(), mFlangeEnd->GetPos().y(), mFlangeEnd->GetPos().z());
    application.AddTypicalCamera(targ + core::vector3df(+0.1f, +0.2f, -0.2f), targ);

    // Initial setup
    application.AssetBindAll();
    application.AssetUpdateAll();
    my_system.SetupInitial();


    // Solver and timestepper
    my_system.SetSolverType(ChSolver::Type::MINRES);
    my_system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(460);
    my_system.SetMaxItersSolverStab(460);
    my_system.SetTolForce(1e-13);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);
#ifdef USE_NSC
    auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
    auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
    my_system.SetStabSolver(ip_solver_stab);
    my_system.SetSolver(ip_solver_speed);
    ip_solver_speed->GetEngine()->SetSparsityPatternLock(true);
    ip_solver_speed->GetEngine()->ForceSparsityPatternUpdate(true);
    // ip_solver_speed->SetVerbose(true);
    application.GetSystem()->Update();
#endif

    // Loop
    application.SetPlotLinkFrames(true);
    auto timestep = contact_model_NSC ? 0.002 : 0.001;

    application.SetTimestep(timestep);
    application.SetVideoframeSaveInterval(1);
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        application.DoStep();

        if (!application.GetPaused()) {
#ifdef USE_NSC
            std::cout << std::setprecision(8)
                << "RotAngle: " << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << "; "
                << std::endl;
#else
            std::cout << "RotAngle: " << rotMot->GetMotorRot()*CH_C_RAD_TO_DEG << std::endl;
#endif
        }
        application.EndScene();
    }

    return 0;
}

