//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010-2012 Alessandro Tasora
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//
//   Demo code about
//
//     - collisions and contacts
//     - sharing a ChMaterialSurface property between bodies
//
//       (This is just a possible method of integration
//       of Chrono::Engine + Irrlicht: many others
//       are possible.)
//
//	 CHRONO
//   ------
//   Multibody dinamics engine
//
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
///////////////////////////////////////////////////

#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/assets/ChTexture.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include <chrono_interiorpoint/ChInteriorPoint.h>

// Use the namespace of Chrono

using namespace chrono;

// Use the main namespaces of Irrlicht
using namespace irr;

using namespace core;
using namespace scene;
using namespace video;
using namespace io;
using namespace gui;

std::shared_ptr<ChBody> generic_body_ptr;


// Create a bunch of ChronoENGINE rigid bodies that
// represent bricks in a large wall.
void create_system(ChSystem& mphysicalSystem) {

	// Create a material that will be shared between bricks
	auto mmaterial = std::make_shared<ChMaterialSurface>();
	mmaterial->SetFriction(0.4f);
	
	if (false) // material selector
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

	// Create the floor using
	// fixed rigid body of 'box' type:
	auto mrigidFloor = std::make_shared<ChBodyEasyBox>(250, 4, 250,  // x,y,z size
		1000,         // density
		true,         // collide enable?
		true);       // visualization?
	mrigidFloor->SetPos(ChVector<>(0, -2, 0));
	mrigidFloor->SetMaterialSurface(mmaterial);
	mrigidFloor->SetBodyFixed(true);

	mphysicalSystem.Add(mrigidFloor);


	
	if (false)
	{
		// Create bricks
		auto mrigidBody = std::make_shared<ChBodyEasyBox>(4, 2, 2,  // x,y,z size
			100,         // density
			true,        // collide enable?
			true);      // visualization?
		mrigidBody->SetPos(ChVector<>(0, 4, 0));
		mrigidBody->SetMaterialSurface(mmaterial);  // use shared surface properties

        generic_body_ptr = mrigidBody;
		mphysicalSystem.Add(mrigidBody);

		// optional, attach a texture for better visualization
		auto mtexture = std::make_shared<ChTexture>();
		mtexture->SetTextureFilename(GetChronoDataFile("cubetexture_borders.png"));
		mrigidBody->AddAsset(mtexture);
	}
	
	if (true)
	{
		// Create a ball that will collide with wall
		auto mrigidBall = std::make_shared<ChBodyEasySphere>(1,       // radius
			8000,    // density
			true,    // collide enable?
			true);  // visualization?
		mrigidBall->SetMaterialSurface(mmaterial);
		mrigidBall->SetPos(ChVector<>(0, 3, 0));
		mrigidBall->SetPos_dt(ChVector<>(0, -1, 0));          // set initial speed

		mphysicalSystem.Add(mrigidBall);
        generic_body_ptr = mrigidBall;

		// optional, attach a texture for better visualization
		auto mtextureball = std::make_shared<ChTexture>();
		mtextureball->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
		mrigidBall->AddAsset(mtextureball);
	}

	
}


int main(int argc, char* argv[]) {
    // Create a ChronoENGINE physical system
    ChSystem mphysicalSystem;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    irrlicht::ChIrrApp application(&mphysicalSystem, L"Balls in bucket", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 290, 190);
    irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-5, 5, -5), core::vector3df(0, 0, 0));
    //ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15, 14, -30), core::vector3df(0, 5, 0));

    //
    // HERE YOU POPULATE THE MECHANICAL SYSTEM OF CHRONO...
    //

    // Create all the rigid bodies.
    create_system(mphysicalSystem);

    // Use this function for adding a ChIrrNodeAsset to all items
    // If you need a finer control on which item really needs a visualization proxy in
    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.
    application.AssetBindAll();

    // Use this function for 'converting' into Irrlicht meshes the assets
    // into Irrlicht-visualizable meshes
    application.AssetUpdateAll();

    // Prepare the physical system for the simulation

    mphysicalSystem.SetSolverType(ChSolver::Type::SOR_MULTITHREAD);

    //mphysicalSystem.SetUseSleeping(true);

    mphysicalSystem.SetMaxPenetrationRecoverySpeed(1.6);  // used by Anitescu stepper only
    mphysicalSystem.SetMaxItersSolverSpeed(40);
    mphysicalSystem.SetMaxItersSolverStab(20);  // unuseful for Anitescu, only Tasora uses this
    mphysicalSystem.SetSolverWarmStarting(true);
    mphysicalSystem.SetParallelThreadNumber(4);

    // Change solver to IP
    auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
    auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
    mphysicalSystem.SetStabSolver(ip_solver_stab);
    mphysicalSystem.SetSolver(ip_solver_speed);
    application.GetSystem()->Update();

    //// Change solver to Matlab external linear solver, for max precision in benchmarks
    //ChMatlabEngine matlab_engine;
    //ChLcpMatlabSolver* matlab_solver_stab = new ChLcpMatlabSolver(matlab_engine);
    //ChLcpMatlabSolver* matlab_solver_speed = new ChLcpMatlabSolver(matlab_engine);
    //mphysicalSystem.SetStabSolver(matlab_solver_stab);
    //mphysicalSystem.SetSolver(matlab_solver_speed);
    //application.GetSystem()->Update();


    //
    // THE SOFT-REAL-TIME CYCLE
    //

    application.SetStepManage(true);
    application.SetTimestep(0.02);

    int step_counter = 0;
    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

        irrlicht::ChIrrTools::drawGrid(application.GetVideoDriver(), 5, 5, 20, 20,
                             ChCoordsys<>(ChVector<>(0, 0.04, 0), Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
                             video::SColor(50, 90, 90, 150), true);

        application.DrawAll();

        application.DoStep();

        step_counter++;
        std::cout << "Step " << step_counter << std::endl;
        //GetLog()
        //    << "Force:" << generic_body_ptr->Get_Xforce()
        //    << "Speed: " << generic_body_ptr->GetPos_dt()
        //    << "ContactForce: " << generic_body_ptr->GetContactForce()
        //    << "\n";


        application.GetVideoDriver()->endScene();
    }

    return 0;
}
