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

#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono_cvxopt/ChSolverCvxoptConeQp.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "solver/ChSolverBB.h"

// Use the namespace of Chrono
using namespace irr;
using namespace irr::core;
using namespace irr::scene;
using namespace irr::video;
using namespace irr::io;
using namespace irr::gui;

using namespace chrono;

int main(int argc, char* argv[]) {
    // Create a ChronoENGINE physical system
    ChSystemNSC mphysicalSystem;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    irrlicht::ChIrrApp application(&mphysicalSystem, L"CVXOPT: friction", core::dimension2d<u32>(800, 600), false,
                                   true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f),
                                              core::vector3df(30.f, 80.f, 60.f), 290, 190);
    irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-5, 5, -5),
                                              core::vector3df(0, 0, 0));
    // ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15, 14, -30), core::vector3df(0, 5,
    // 0));

    const int num_bodies = 5;
    const double body_mass = 10;
    const double fric_coeff = 0.5;
    const double gravity = 10;

    // const double force_max = 0;
    // const double incl_max = CH_C_PI / 3;

    const double force_max = body_mass * gravity * fric_coeff * 1.5;
    const double incl_max = 0;

    std::array<double, num_bodies> forces_additional;
    std::array<double, num_bodies> inclination;
    std::cout << std::endl;
    std::cout << "Fric coeff: " << fric_coeff << " (angle: " << atan(fric_coeff) * 180.0 / CH_C_PI << " deg)"
              << std::endl;
    std::cout << "Body weight: " << body_mass * gravity << " N" << std::endl;
    std::cout << "Expected max friction: " << body_mass * gravity * fric_coeff << " N" << std::endl;

    for (auto body_sel = 0; body_sel < num_bodies; ++body_sel) {
        inclination[body_sel] = body_sel * incl_max / std::max(num_bodies - 1, 1);
        forces_additional[body_sel] = body_sel * force_max / std::max(num_bodies - 1, 1);
        std::cout << "Body " << body_sel << " | incl: " << inclination[body_sel] * 180.0 / CH_C_PI
                  << " deg; force horiz.: " << forces_additional[body_sel] << " N" << std::endl;
    }

    application.GetSystem()->Set_G_acc(ChVector<>(0, -gravity, 0));

    // Material
    auto mmaterial = std::make_shared<ChMaterialSurfaceNSC>();
    mmaterial->SetCompliance(0.0f);
    mmaterial->SetComplianceT(0.0f);
    mmaterial->SetDampingF(0.0f);
    mmaterial->SetFriction(fric_coeff);

    // Asset
    auto floorColor = std::make_shared<ChColorAsset>();
    floorColor->SetColor(ChColor(0.25, 0.25, 0.25));
    const double box_edge = 1.0;

    for (auto body_sel = 0; body_sel < num_bodies; ++body_sel) {
        // create the floor
        auto mrigidFloor = std::make_shared<ChBodyEasyBox>(2, 1, 10, 1000, true, true, ChMaterialSurface::NSC);
        mrigidFloor->SetPos(ChMatrix33<>(Q_from_AngAxis(inclination[body_sel], VECT_X)) * ChVector<>(0, -0.5, 0) +
                            ChVector<>(body_sel * 2.0, 0, 0));
        mrigidFloor->SetRot(Q_from_AngAxis(inclination[body_sel], VECT_X));
        mrigidFloor->SetMaterialSurface(mmaterial);
        mrigidFloor->SetBodyFixed(true);
        mrigidFloor->AddAsset(floorColor);
        mphysicalSystem.Add(mrigidFloor);

        // give color
        auto color = std::make_shared<ChColorAsset>();
        color->SetColor(ChColor::ComputeRainbowColor(body_sel, num_bodies));

        // Create a box that will collides with the floor
        auto mrigidBox = std::make_shared<ChBodyEasyBox>(1, 1, 1, 1000, true, true, ChMaterialSurface::NSC);
        mrigidBox->SetRot(Q_from_AngAxis(inclination[body_sel], VECT_X));
        mrigidBox->SetPos(ChMatrix33<>(Q_from_AngAxis(inclination[body_sel], VECT_X)) *
                              ChVector<>(0, +0.5 * box_edge, 0) +
                          ChVector<>(body_sel * box_edge * 2, 0, 0));
        mrigidBox->Set_Scr_force(VECT_Z * forces_additional[body_sel]);
        mrigidBox->Set_Scr_torque(mrigidBox->Get_Scr_force().Cross(VECT_Y * 0.5 * box_edge));
        mrigidBox->SetMaterialSurface(mmaterial);
        mrigidBox->SetMass(body_mass);
        mrigidBox->AddAsset(color);

        mphysicalSystem.Add(mrigidBox);

        //// Apply forces to body
        // auto force = std::make_shared<ChForce>();
        // mrigidBox->AddForce(force);
        // force->SetMode(ChForce::FORCE);
        // force->SetFrame(ChForce::WORLD);
        // force->SetDir(VECT_Z);
        // force->SetMforce(forces_additional[body_sel]);
    }

    application.AssetBindAll();
    application.AssetUpdateAll();

    // Set CVXOPT ConeQP solver
    auto ip_solver_stab = std::make_shared<ChSolverCvxoptConeQp>();
    auto ip_solver_speed = std::make_shared<ChSolverCvxoptConeQp>();
    mphysicalSystem.SetStabSolver(ip_solver_stab);
    mphysicalSystem.SetSolver(ip_solver_speed);
    ip_solver_speed->SetSparsityPatternLock(true);
    ip_solver_speed->ForceSparsityPatternUpdate(true);
    // ip_solver_speed->SetVerbose(true);
    application.GetSystem()->Update();

    //// Set CVXOPT ConeQP solver
    // auto solver_stab = std::make_shared<ChSolverBB>();
    // auto solver_speed = std::make_shared<ChSolverBB>();
    // mphysicalSystem.SetStabSolver(solver_stab);
    // mphysicalSystem.SetSolver(solver_speed);
    ////ip_solver_speed->SetVerbose(true);
    // application.GetSystem()->Update();

    application.SetStepManage(true);
    application.SetTimestep(0.01);
    application.SetTryRealtime(true);
    // application.SetPaused(true);
    application.SetContactsDrawMode(irrlicht::ChIrrTools::eCh_ContactsDrawMode::CONTACT_FORCES);
    application.SetContactsLabelMode(irrlicht::ChIrrTools::eCh_ContactsLabelMode::CONTACT_FORCES_N_VAL);
    int step_counter = 0;
    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

        application.DrawAll();

        if (!application.GetPaused()) {
            application.DoStep();
            std::cout << "Gap: " << std::static_pointer_cast<ChSolverCvxoptConeQp>(application.GetSystem()->GetSolver())->GetEngine().GetGap() << "; "
                    << "RelGap: " << std::static_pointer_cast<ChSolverCvxoptConeQp>(application.GetSystem()->GetSolver())->GetEngine().GetRelativeGap() << "; "
                    << "PrimInfeas: " << std::static_pointer_cast<ChSolverCvxoptConeQp>(application.GetSystem()->GetSolver())->GetEngine().GetPrimalInfeasibility() << "; "
                    << "DualInfeas: " << std::static_pointer_cast<ChSolverCvxoptConeQp>(application.GetSystem()->GetSolver())->GetEngine().GetDualInfeasibility() << "; "
                    << std::endl;
        }


        application.GetVideoDriver()->endScene();
    }

    return 0;
}