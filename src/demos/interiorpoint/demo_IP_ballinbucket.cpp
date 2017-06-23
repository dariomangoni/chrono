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

// Create a bunch of ChronoENGINE rigid bodies that
// represent bricks in a large wall.

void create_bucket(ChSystemNSC& mphysicalSystem) {

    // Create a material that will be shared between bricks
    auto mmaterial = std::make_shared<ChMaterialSurfaceNSC>();

    mmaterial->SetFriction(0.4f);
    mmaterial->SetCompliance(0.0000005f);
    mmaterial->SetComplianceT(0.0000005f);
    mmaterial->SetDampingF(0.2f);

    // Create bricks
	int box_edges = 36;
	double wall_thick = 0.01;
	double box_apothem = 0.5;
	double box_height = 0.25;
	double ball_radius = 0.035;
	double ball_clearance = 0.005;
    int ball_arrays = 5;


	
	double alfa = 2*PI/box_edges;
	double width = 2*box_apothem*tan(alfa/2); // width bricks

	// Create Polyhedron
	for (int edge_k = 0; edge_k < box_edges; edge_k++)
	{
		double alfa_k = alfa*edge_k;
        auto wall = std::make_shared<ChBodyEasyBox>(width, box_height, wall_thick,  // x,y,z size
                                                    1000,         // density
                                                    true,         // collide enable?
                                                    (alfa_k >= 0 && alfa_k < PI) ? true : false);

		ChQuaternion<double> quat(cos(alfa_k/2), 0, sin(alfa_k / 2), 0);
		wall->SetRot(quat);
		wall->SetPos(ChVector<>((box_apothem + wall_thick / 2)*sin(alfa_k), box_height / 2, (box_apothem + wall_thick / 2)*cos(alfa_k)));
		wall->SetMaterialSurface(mmaterial);
		wall->SetBodyFixed(true);
		mphysicalSystem.Add(wall);

	}


    double square_box_edge = box_apothem*2;
    auto mtextureball = std::make_shared<ChTexture>();
    mtextureball->SetTextureFilename(GetChronoDataFile("bluwhite.png"));
    for (auto ball_array = 0; ball_array < ball_arrays; ++ball_array)
        for (auto body_pos_x = -square_box_edge / 2 + ball_clearance; body_pos_x <= square_box_edge / 2 - ball_clearance; body_pos_x += 2 * (ball_radius + ball_clearance))
            for (auto body_pos_z = -square_box_edge / 2 + ball_clearance; body_pos_z <= square_box_edge / 2 - ball_clearance; body_pos_z += 2 * (ball_radius + ball_clearance))
            {
                if (sqrt(body_pos_x*body_pos_x + body_pos_z*body_pos_z) > box_apothem - ball_clearance)
                    continue;

                auto mrigidBall = std::make_shared<ChBodyEasySphere>(ball_radius, 8000, true, true);  // visualization?
                mrigidBall->SetPos(ChVector<>(body_pos_x, (2 * ball_array + 1)*(ball_radius + ball_clearance), body_pos_z));
                mrigidBall->AddAsset(mtextureball);
                mrigidBall->GetMaterialSurfaceNSC()->SetFriction(0.4f);  // use own (not shared) matrial properties
                mrigidBall->GetMaterialSurfaceNSC()->SetCompliance(0.0);
                mrigidBall->GetMaterialSurfaceNSC()->SetComplianceT(0.0);
                mrigidBall->GetMaterialSurfaceNSC()->SetDampingF(0.2f);
                mrigidBall->AddAsset(mtextureball);
                mphysicalSystem.Add(mrigidBall);
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

    
}


int main(int argc, char* argv[]) {
    // Create a ChronoENGINE physical system
    ChSystemNSC mphysicalSystem;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
	irrlicht::ChIrrApp application(&mphysicalSystem, L"Balls in bucket", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
	irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
	irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice(), core::vector3df(70.f, 120.f, -90.f), core::vector3df(30.f, 80.f, 60.f), 290, 190);
	irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-0.5, 0.5, -0.5), core::vector3df(0, 0, 0));
    //ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(-15, 14, -30), core::vector3df(0, 5, 0));

    //
    // HERE YOU POPULATE THE MECHANICAL SYSTEM OF CHRONO...
    //

    // Create all the rigid bodies.
    create_bucket(mphysicalSystem);

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
    ip_solver_speed->SetNullPivotDetection(true, 1e-18);
	ip_solver_speed->SetUseSymmetry(true);
    ip_solver_speed->SetVerbose(true);

	application.GetSystem()->Update();




    //
    // THE SOFT-REAL-TIME CYCLE
    //

    application.SetStepManage(true);
    application.SetTimestep(0.02);
	//application.SetPaused(true);

	ChTimer<> timer;
	timer.start();
    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

		irrlicht::ChIrrTools::drawGrid(application.GetVideoDriver(), 5, 5, 20, 20,
                             ChCoordsys<>(ChVector<>(0, 0.04, 0), Q_from_AngAxis(CH_C_PI / 2, VECT_X)),
                             video::SColor(50, 90, 90, 150), true);

        application.DrawAll();

        application.DoStep();

        application.GetVideoDriver()->endScene();

        if (ip_solver_speed->GetSolverCalls() > 25)
            break;

    }
	timer.stop();


    std::cout << "Time spent: " << ip_solver_speed->GetTimeSolve_Assembly() + ip_solver_speed->GetTimeSolve_SolverCall()  << "; IP calls: " << ip_solver_speed->GetIPSolverCalls() << "; IP iterations: " << ip_solver_speed->GetIPIterations() << std::endl;
    std::cout << "Mass matrix size: " << ip_solver_speed->GetMassMatrixDimension() << "; Inequalities: " << ip_solver_speed->GetUnilateralConstraintsMatrixRows()<< "; Equalities: " << ip_solver_speed->GetBilateralConstraintsMatrixRows() <<std::endl;
	std::cout << "Time spent (external timer): " << timer.GetTimeSeconds() << std::endl;

	getchar();

    return 0;
}
