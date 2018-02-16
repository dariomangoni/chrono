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


void run_test(int total_balls)
{
	double ball_radius = 0.5;
	double min_ball_wall_clearance = 0.001;

	// Create a ChronoENGINE physical system
	ChSystemNSC mphysicalSystem;

	// Create the Irrlicht visualization (open the Irrlicht device,
	// bind a simple user interface, etc. etc.)
	irrlicht::ChIrrApp application(&mphysicalSystem, L"Critical Ball", core::dimension2d<u32>(1080, 1080), false, true);

	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
	irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
	irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, -120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 290, 190);
	irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(+1.0, +0.5, +1.0), core::vector3df(0, 0, 0));
	//ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15, 14, -30), core::vector3df(0, 5, 0));

	//
	// HERE YOU POPULATE THE MECHANICAL SYSTEM OF CHRONO...
	//

	// Create a material that will be shared between bricks
	auto hard_material = std::make_shared<ChMaterialSurfaceNSC>();
	hard_material->SetFriction(0.4f);
	hard_material->SetCompliance(0.0f);
	hard_material->SetComplianceT(0.0f);
	hard_material->SetDampingF(0.0f);

	// Create the floor using
	// fixed rigid body of 'box' type:
	auto mrigidFloor = std::make_shared<ChBodyEasyBox>(250, 1, 250, 1000, true, true);
	mrigidFloor->SetPos(ChVector<>(0, -0.5, 0));
	mrigidFloor->SetMaterialSurface(hard_material);
	mrigidFloor->SetBodyFixed(true);
	mphysicalSystem.Add(mrigidFloor);

	auto mtextureball = std::make_shared<ChTexture>();
	mtextureball->SetTextureFilename(GetChronoDataFile("bluwhite.png"));

	// Create a ball that will collide with wall
	auto fixed_truss = std::make_shared<ChBodyEasyBox>(0.1, 0.1, 0.1, 1000, false, false, ChMaterialSurface::NSC);
	fixed_truss->SetBodyFixed(true);
	fixed_truss->SetPos(VNULL);
	mphysicalSystem.Add(fixed_truss);

	double density = 7500;
	double clearance = 0.01;
	auto ball_rad = [density](double mass) { return std::pow(mass* 3.0 / 4.0 / density / CH_C_PI, 1.0 / 3.0); };
	auto mass_choice = [density](int index) { return std::pow(10.0, static_cast<double>(index)); };
	//auto mass_choice = [](int) { return 1e12; };

	double height = 0;
	for (auto ball_sel = 0; ball_sel < total_balls; ++ball_sel)
	{
		height += clearance + ball_rad(mass_choice(ball_sel));
		auto ball_temp = std::make_shared<ChBodyEasySphere>(ball_rad(mass_choice(ball_sel)), density, true, true);
		ball_temp->SetMaterialSurface(hard_material);
		ball_temp->SetPos(ChVector<>(0, height, 0));
		mphysicalSystem.Add(ball_temp);
		ball_temp->AddAsset(mtextureball);

		//std::cout << "Height: " << height << std::endl;

		//std::cout << "Ball" << ball_sel << ": " << ball_temp->GetMass() << " kg" << std::endl;

		auto ball_temp_link = std::make_shared<ChLinkLockPrismatic>();
		ball_temp_link->Initialize(ball_temp, fixed_truss, ChCoordsys<>(ChVector<>(0.0, -1.0, 0.0), Q_from_AngX(-CH_C_PI_2)));
		mphysicalSystem.AddLink(ball_temp_link);

		height += ball_rad(mass_choice(ball_sel));


	}


	application.AssetBindAll();
	application.AssetUpdateAll();


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
	ip_solver_speed->RecordHistory(false);
	//ip_solver_speed->SetVerbose(true);
	application.GetSystem()->Update();


	//
	// THE SOFT-REAL-TIME CYCLE
	//

	application.SetStepManage(true);
	application.SetTimestep(0.05);
	//application.SetTryRealtime(true);
	//application.SetPaused(true);

	int step_counter = 0;
	while (application.GetDevice()->run()) {
		application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

		irrlicht::ChIrrTools::drawGrid(application.GetVideoDriver(), 5, 5, 20, 20,
			ChCoordsys<>(ChVector<>(0, 0.04, 0), Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
			video::SColor(50, 90, 90, 150), true);

		application.DrawAll();

		application.DoStep();



		application.GetVideoDriver()->endScene();

		if (ip_solver_speed->GetSolverCalls() > 50)
			break;
	}

	std::cout << "Time spent: " << ip_solver_speed->GetTimeSolve_Assembly() + ip_solver_speed->GetTimeSolve_SolverCall() << "; IP calls: " << ip_solver_speed->GetIPSolverCalls() << "; IP iterations: " << ip_solver_speed->GetIPIterations() << std::endl;
	std::cout << "Mass matrix size: " << ip_solver_speed->GetMassMatrixDimension() << "; Inequalities: " << ip_solver_speed->GetUnilateralConstraintsMatrixRows() << "; Equalities: " << ip_solver_speed->GetBilateralConstraintsMatrixRows() << std::endl;

}

int main(int argc, char* argv[]) {
	

	//for (auto ball_count = 1; ball_count < 25; ++ball_count)
	//	run_test(ball_count);

	run_test(12);

	getchar();

    return 0;
}
