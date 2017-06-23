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



int main(int argc, char* argv[]) {
	
	double ball_radius = 0.5;
	double min_ball_wall_clearance = 0.001;
	double ball_wall_initial_clearance = 0.1;
	
    // Create a ChronoENGINE physical system
    ChSystemNSC mphysicalSystem;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    irrlicht::ChIrrApp application(&mphysicalSystem, L"Critical Ball", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 290, 190);
    irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(0.0, 5.0, 0.0), core::vector3df(0.0, 0.0, 0.0));
    //ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15, 14, -30), core::vector3df(0, 5, 0));
	application.SetPlotLinkFrames(true);

	//// build the orthographic projection matrix
	//// play with the values here to change the coordinate system
	//// the values I use below, set a coordinate system from (-0.5,-0.5) to (0.5,0.5). This means that the center of the screen is at (0,0)
	//irr::core::matrix4 orthoProjection;
	//orthoProjection.buildProjectionMatrixOrthoLH(1.0f, 1.0f, 0.0f, 1.0f);

	//application.GetSceneManager()->getActiveCamera()->setProjectionMatrix(orthoProjection, true);

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
	auto mrigidFloor = std::make_shared<ChBodyEasyBox>(250, 1, 250, 1000, true,	true, ChMaterialSurface::NSC);
	mrigidFloor->SetPos(ChVector<>(0, -0.5, 0));
	mrigidFloor->SetMaterialSurface(hard_material);
	mrigidFloor->SetBodyFixed(true);

	mphysicalSystem.Add(mrigidFloor);

	// Create a ball that will collide with wall
	auto mrigidBall = std::make_shared<ChBodyEasySphere>(ball_radius, 8000, true, false, ChMaterialSurface::NSC);
	mrigidBall->SetMaterialSurface(hard_material);
	mrigidBall->SetMass(10);
	mrigidBall->SetPos(ChVector<>(0, ball_radius*1.01, 0));
	mphysicalSystem.Add(mrigidBall);

	auto mtextureball = std::make_shared<ChTexture>();
	mtextureball->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
	mrigidBall->AddAsset(mtextureball);


	// Create bricks
	const int box_edges = 3;
	double wall_thick = 0.01;
	double box_apothem = ball_radius + ball_wall_initial_clearance;
	double box_height = ball_radius;

	  
	auto fixed_truss = std::make_shared<ChBodyEasyBox>(0.1, 0.1, 0.1, 1000, false, false, ChMaterialSurface::NSC);
	fixed_truss->SetBodyFixed(true);
	fixed_truss->SetPos(VNULL);
	mphysicalSystem.Add(fixed_truss);

	double alfa = 2 * PI / box_edges;
	double width = 2 * ball_radius*0.99*tan(alfa / 2); // width bricks

	std::cout << alfa * 180 / CH_C_PI << std::endl;

	ChLinkLimit wall_link_lim[box_edges];

	// Create Polyhedron
	for (int edge_k = 0; edge_k < box_edges; edge_k++)
	{

		double alfa_k = alfa*edge_k;
		auto wall = std::make_shared<ChBodyEasyBox>(width, box_height, wall_thick,  // x,y,z size
			1000,         // density
			true,         // collide enable?
			true);

		ChQuaternion<double> quat(cos(alfa_k / 2), 0, sin(alfa_k / 2), 0);
		wall->SetRot(quat);
		wall->SetPos(ChVector<>((box_apothem + wall_thick / 2)*sin(alfa_k), box_height / 2, (box_apothem + wall_thick / 2)*cos(alfa_k)));
		wall->SetMaterialSurface(hard_material);
		mphysicalSystem.Add(wall);

		auto wall_link = std::make_shared<ChLinkLockLock>();
		wall_link->Initialize(wall, fixed_truss, ChCoordsys<>(VNULL, Q_from_AngY(alfa_k)));
		auto wall_link_function = std::make_shared<ChFunction_Ramp>();
		wall_link_function->Set_ang(+1.0);
		wall_link_function->Set_y0(0.0);
		wall_link_lim[edge_k].Set_min(-0.5);
		wall_link_lim[edge_k].Set_max(+0.5);
		wall_link->SetLimit_Z(&wall_link_lim[edge_k]);
		wall_link->SetMotion_Z(wall_link_function);
		mphysicalSystem.AddLink(wall_link);

	}



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
    ip_solver_speed->RecordHistory(false);
    ip_solver_speed->SetVerbose(true);
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
    //application.SetTryRealtime(true);
	application.SetPaused(true);

    int step_counter = 0;
    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

        irrlicht::ChIrrTools::drawGrid(application.GetVideoDriver(), 5, 5, 20, 20,
                             ChCoordsys<>(ChVector<>(0, 0.04, 0), Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
                             video::SColor(50, 90, 90, 150), true);

        application.DrawAll();

        application.DoStep();


        application.GetVideoDriver()->endScene();

        //if (ip_solver_speed->GetSolverCalls() > 100)
        //    break;
    }

	std::cout << "Time spent: " << ip_solver_speed->GetTimeSolve_Assembly() + ip_solver_speed->GetTimeSolve_SolverCall() << "; IP calls: " << ip_solver_speed->GetIPSolverCalls() << "; IP iterations: " << ip_solver_speed->GetIPIterations() << std::endl;
	std::cout << "Mass matrix size: " << ip_solver_speed->GetMassMatrixDimension() << "; Inequalities: " << ip_solver_speed->GetUnilateralConstraintsMatrixRows() << "; Equalities: " << ip_solver_speed->GetBilateralConstraintsMatrixRows() << std::endl;
	getchar();

    return 0;
}
