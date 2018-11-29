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
#include "chrono/assets/ChTexture.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_cvxopt/ChSolverCvxoptConeQp.h"

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
    irrlicht::ChIrrApp application(&mphysicalSystem, L"CVXOPT: ball", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 290, 190);
    irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-5, 5, -5), core::vector3df(0, 0, 0));
    application.SetContactsDrawMode(irrlicht::ChIrrTools::eCh_ContactsDrawMode::CONTACT_FORCES);
    application.SetContactsLabelMode(irrlicht::ChIrrTools::eCh_ContactsLabelMode::CONTACT_FORCES_N_VAL);
    application.SetSymbolscale(1e-2);

    // Create a material that will be shared between bricks
    auto mmaterial = std::make_shared<ChMaterialSurfaceNSC>();
    mmaterial->SetFriction(1.0f);

    if (true) // material selector
    {
        mmaterial->SetCompliance(0.0000005f);
        mmaterial->SetComplianceT(0.0000005f);
        mmaterial->SetDampingF(0.5f);
    }
    else
    {
        mmaterial->SetCompliance(0.0f);
        mmaterial->SetComplianceT(0.0f);
        mmaterial->SetDampingF(0.0f);
    }

    // Create the floor using 'box':
    auto mrigidFloor = std::make_shared<ChBodyEasyBox>(250, 4, 250, 1000, true, true);
    mrigidFloor->SetPos(ChVector<>(0, -2, 0));
    mrigidFloor->SetMaterialSurface(mmaterial);
    mrigidFloor->SetBodyFixed(true);

    mphysicalSystem.Add(mrigidFloor);

    // Ball with collision (to test unilateral (cone) constraints)
    auto mrigidBall = std::make_shared<ChBodyEasySphere>(1, 8000, true, true);  // visualization?
    mrigidBall->SetMaterialSurface(mmaterial);
    mrigidBall->SetMass(10);
    mrigidBall->SetPos(ChVector<>(0, 1.5, 0));

    auto mtextureball = std::make_shared<ChTexture>();
    mtextureball->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    mrigidBall->AddAsset(mtextureball);

    mphysicalSystem.Add(mrigidBall);






    // Box, no collision (to test bilateral constraints)
    auto mrigidBox = std::make_shared<ChBodyEasyBox>(1, 1, 1, 1000, false, true, ChMaterialSurface::NSC);
    mrigidBox->SetMaterialSurface(mmaterial);
    mrigidBox->SetMass(10);
    mrigidBox->SetPos(ChVector<>(1.0, 5, 0));

    auto color = std::make_shared<ChColorAsset>();
    color->SetColor(ChColor(0.0, 0.0, 1.0));
    mrigidBox->AddAsset(color);

    mphysicalSystem.Add(mrigidBox);

    auto link = std::make_shared<ChLinkLockRevolute>();
    link->Initialize(mrigidFloor, mrigidBox, ChCoordsys<>(ChVector<>(0.0,5,0.0), QUNIT));
    mphysicalSystem.Add(link);



    application.AssetBindAll();
    application.AssetUpdateAll();

    // Change solver to IP
    auto ip_solver_stab = std::make_shared<ChSolverCvxoptConeQp>();
    auto ip_solver_speed = std::make_shared<ChSolverCvxoptConeQp>();
    mphysicalSystem.SetStabSolver(ip_solver_stab);
    mphysicalSystem.SetSolver(ip_solver_speed);
    ip_solver_speed->SetSparsityPatternLock(true);
    ip_solver_speed->ForceSparsityPatternUpdate(true);
    //ip_solver_speed->SetVerbose(true);
    application.GetSystem()->Update();


    application.GetSystem()->SetMaxPenetrationRecoverySpeed(1e-2);

    application.SetStepManage(true);
    application.SetTimestep(0.01);
    //application.SetTryRealtime(true);
    //application.SetPaused(true);

    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

        irrlicht::ChIrrTools::drawGrid(application.GetVideoDriver(), 5, 5, 20, 20,
            ChCoordsys<>(ChVector<>(0, 0.04, 0), Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
            video::SColor(50, 90, 90, 150), true);

        application.DrawAll();

        application.DoStep();

        application.GetVideoDriver()->endScene();
    }

    getchar();

    return 0;
}