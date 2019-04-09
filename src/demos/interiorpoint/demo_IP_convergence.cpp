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
#include "chrono/core/ChRealtimeStep.h"
#include "chrono/solver/ChSolverMINRES.h"

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

// Static data used for this simple demo

std::vector<std::shared_ptr<ChBody>> mspheres;

float STATIC_COMPLIANCE = 0.1f * (10.0f / 1000) / 5;  // as 1/K, in m/N. es: 10mm/500N

void create_items(irrlicht::ChIrrAppInterface& application, std::shared_ptr<ChMaterialSurface> material) {
    // Create some spheres in a vertical stack

    bool do_wall = false;
    bool do_stack = true;
    bool do_oddmass = true;
    bool do_spheres = true;
    bool do_heavyonside = true;

    double sphrad = 0.2;
    double dens = 1000;
    double sphmass = dens * (4. / 3.) * CH_C_PI * pow(sphrad, 3);
    double sphinertia = (2. / 5.) * sphmass * pow(sphrad, 2);

    if (do_stack) {
        int nbodies = 15;

        double totmass = 0;
        double level = 0;
        double sphrad_base = 0.2;
        double oddfactor = 100;

        for (int bi = 0; bi < nbodies; bi++)  // N. of vert. bricks
        {
            double sphrad = sphrad_base;
            if (do_oddmass && bi == (nbodies - 1))
                sphrad = sphrad * pow(oddfactor, 1. / 3.);
            double dens = 1000;

            std::shared_ptr<ChBody> mrigidBody;

            if (do_spheres) {
                mrigidBody = std::make_shared<ChBodyEasySphere>(sphrad,  // radius
                                                                dens,    // density
                                                                true,    // collide enable?
                                                                true);  // visualization?

                mrigidBody->SetPos(ChVector<>(0.5, sphrad + level, 0.7));
                mrigidBody->AddAsset(std::make_shared<ChTexture>(GetChronoDataFile("bluwhite.png")));

                application.GetSystem()->Add(mrigidBody);
            } else {
                mrigidBody = std::make_shared<ChBodyEasyBox>(sphrad, sphrad, sphrad,  // x,y,z size
                                                                          dens,                    // density
                                                                          true,                    // collide enable?
                                                                          true);                  // visualization?
                mrigidBody->SetPos(ChVector<>(0.5, sphrad + level, 0.7));
                mrigidBody->AddAsset(std::make_shared<ChTexture>(GetChronoDataFile("cubetexture_bluwhite.png")));

                application.GetSystem()->Add(mrigidBody);
            }

			mrigidBody->SetMaterialSurface(material);


            mspheres.push_back(mrigidBody);

            level += sphrad * 2;
            totmass += mrigidBody->GetMass();
        }

        GetLog() << "Expected contact force at bottom F=" << (totmass * application.GetSystem()->Get_G_acc().y()) << "\n";
    }

    if (do_wall)
        for (int ai = 0; ai < 1; ai++)  // N. of walls
        {
            for (int bi = 0; bi < 10; bi++)  // N. of vert. bricks
            {
                for (int ui = 0; ui < 15; ui++)  // N. of hor. bricks
                {
                    auto mrigidWall = std::make_shared<ChBodyEasyBox>(3.96, 2, 4,  // radius
                                                                            dens,        // density
                                                                            true,        // collide enable?
                                                                            true);      // visualization?
                    mrigidWall->SetPos(ChVector<>(-8 + ui * 4.0 + 2 * (bi % 2), 1.0 + bi * 2.0, -5 + ai * 6));
					
					mrigidWall->SetMaterialSurface(material);

                    mrigidWall->AddAsset(
                        std::make_shared<ChTexture>(GetChronoDataFile("cubetexture_bluwhite.png")));

                    application.GetSystem()->Add(mrigidWall);
                }
            }
        }

    if (do_heavyonside) {
        double sphrad = 0.2;
        double dens = 1000;
        double hfactor = 100;

        auto mrigidHeavy = std::make_shared<ChBodyEasySphere>(sphrad,          // radius
                                                                       dens * hfactor,  // density
                                                                       true,            // collide enable?
                                                                       true);          // visualization?
        mrigidHeavy->SetPos(ChVector<>(0.5, sphrad + 0.1, -1));
        mrigidHeavy->AddAsset(std::make_shared<ChTexture>(GetChronoDataFile("pinkwhite.png")));

		mrigidHeavy->SetMaterialSurface(material);

        application.GetSystem()->Add(mrigidHeavy);

        GetLog() << "Expected contact deformation at side sphere="
                 << (mrigidHeavy->GetMass() * application.GetSystem()->Get_G_acc().y()) * STATIC_COMPLIANCE << "\n";
    }

    // Create the floor using a fixed rigid body of 'box' type:

    auto mrigidFloor = std::make_shared<ChBodyEasyBox>(50, 4, 50,  // radius
                                                             dens,       // density
                                                             true,       // collide enable?
                                                             true);     // visualization?
    mrigidFloor->SetPos(ChVector<>(0, -2, 0));
    mrigidFloor->SetBodyFixed(true);

	mrigidFloor->SetMaterialSurface(material);

    mrigidFloor->AddAsset(std::make_shared<ChTexture>(GetChronoDataFile("concrete.jpg")));


    application.GetSystem()->Add(mrigidFloor);
}

// Function that forces all spheres in the 'parent' level to be on the same vertical
// axis, without needing any constraint (for simplifying the solver benchmark).
// Also impose no rotation.

void align_spheres(irrlicht::ChIrrAppInterface& application) {
    for (unsigned int i = 0; i < mspheres.size(); ++i) {
        std::shared_ptr<ChBody> body = mspheres[i];
        ChVector<> mpos = body->GetPos();
        mpos.x() = 0.5;
        mpos.z() = 0.7;
        body->SetPos(mpos);
        body->SetRot(QUNIT);
    }
}

int main(int argc, char* argv[]) {
    // Create a ChronoENGINE physical system
    ChSystemNSC mphysicalSystem;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    irrlicht::ChIrrApp application(&mphysicalSystem, L"Critical cases for solver convergence", core::dimension2d<u32>(800, 600),
                         false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    irrlicht::ChIrrWizard::add_typical_Logo(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Sky(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Lights(application.GetDevice());
    irrlicht::ChIrrWizard::add_typical_Camera(application.GetDevice(), core::vector3df(0, 1.5, -3));

	// Create a material that will be shared
	auto material = std::make_shared<ChMaterialSurfaceNSC>();
	material->SetFriction(0.5f);
	//material->SetRestitution(0.0f);
	material->SetCompliance(STATIC_COMPLIANCE);
	material->SetComplianceT(STATIC_COMPLIANCE);
	material->SetDampingF(0.3f);

	// Create all the rigid bodies.
    create_items(application, material);

	/*
    class MyContactCallback : public ChSystemNSC::ChCustomCollisionPointCallback {
      public:
        virtual void ContactCallback(
            const collision::ChCollisionInfo& mcontactinfo,  ///< get info about contact (cannot change it)
            ChMaterialCouple& material)                      ///< you can modify this!
        {
            // Set compliance (normal and tangential at once)
            material.compliance = STATIC_COMPLIANCE;
            material.complianceT = material.compliance;
            material.dampingf = 0.2f;
        };
        ChSystemNSC* msystem;
    };

    MyContactCallback mycontact_callback;           // create the callback object
    mycontact_callback.msystem = &mphysicalSystem;  // will be used by callback
    // Tell the system to use the callback above, per each created contact!
    mphysicalSystem.SetCustomCollisionPointCallback(&mycontact_callback);
	*/

    // Use this function for adding a ChIrrNodeAsset to all already created items (ex. a floor, a wall, etc.)
    // Otherwise use application.AssetBind(myitem); on a per-item basis.
    application.AssetBindAll();
    application.AssetUpdateAll();

    // Modify some setting of the physical system for the simulation, if you want

    mphysicalSystem.SetSolverType(ChSolver::Type::BARZILAIBORWEIN);
    // mphysicalSystem.SetSolverType(ChSystemNSC::LCP_ITERATIVE_SOR);
    mphysicalSystem.SetMaxItersSolverSpeed(60);
    mphysicalSystem.SetMaxItersSolverStab(5);
    mphysicalSystem.SetParallelThreadNumber(1);
	
	// Change solver to Interior-Point
    auto ip_solver_stab = std::make_shared<ChInteriorPoint>();
    auto ip_solver_speed = std::make_shared<ChInteriorPoint>();
    mphysicalSystem.SetStabSolver(ip_solver_stab);
	mphysicalSystem.SetSolver(ip_solver_speed);
	//ip_solver_stab->SetMaxIterations(10);
	//ip_solver_speed->SetMaxIterations(10);
	application.GetSystem()->Update();
	

    // mphysicalSystem.SetUseSleeping(true);

    application.SetStepManage(true);
    application.SetTimestep(0.01);
    application.SetPaused(true);

    //
    // THE SOFT-REAL-TIME CYCLE
    //

    while (application.GetDevice()->run()) {
        application.GetVideoDriver()->beginScene(true, true, SColor(255, 140, 161, 192));

        application.DrawAll();

        align_spheres(application);  // just to simplify test, on y axis only

        application.DoStep();

        application.GetVideoDriver()->endScene();
    }

    return 0;
}
