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
// FEA for shells of Reissner 6-field type
//
// =============================================================================

#include <vector>

#include "chrono/core/ChFileutils.h"

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/timestepper/ChTimestepper.h"

#include "chrono_fea/ChElementShellReissner4.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_fea/ChRotUtils.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_mkl/ChSolverMKL.h"
#include "chrono_postprocess/ChGnuPlot.h"
#include "chrono_fea/ChMeshFileLoader.h"
#include <set>

// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;
using namespace chrono::postprocess;
using namespace irr;

// Output directory

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Create a Chrono::Engine physical system
    ChSystemNSC my_system;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Shells FEA", core::dimension2d<u32>(1024, 768), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.5)), core::vector3dfCH(ChVector<>(0.0, 0.0, 0.0)));
    // application.SetContactsDrawMode(irr::ChIrrTools::CONTACT_DISTANCES);

    application.AddLightWithShadow(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.5)), core::vector3df(0.0, 0.0, 0.0), 0.1, 0.1, 0.5,
        0, 512, video::SColorf((f32)0.8, (f32)0.8, (f32)1.0));


    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // my_system.Set_G_acc(VNULL); or
    my_mesh->SetAutomaticGravity(false);

    // Create a material
    double rho = 7850;
    double E = 200e9;
    double nu = 0.7;

    auto shell_material = std::make_shared<ChMaterialShellReissnerIsothropic>(rho, E, nu);

    std::map<unsigned, std::tuple<std::string, std::vector<unsigned>>> elements_map;
    std::map<unsigned, std::vector<double>> nodes_map;
    std::map<std::string, std::vector<unsigned int>> nset_map;
    std::map<std::string, std::vector<unsigned int>> elset_map;
    std::map<unsigned int, std::shared_ptr<ChNodeFEAxyzrot>> inserted_nodes;

    try {
        ChMeshFileLoader::FromAbaqusFileMOD("C:/workspace/chrono_worktree/data/fea/mesh_good_but_not_perfect.INP", elements_map, nodes_map, nset_map, elset_map);
    }
    catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    auto shell_thickness = 0.01;
    for (auto el_it = elements_map.begin(); el_it!=elements_map.end(); ++el_it)
    {
        if (std::get<0>(el_it->second) == "CPS4")
        {
            auto new_elem = std::make_shared<ChElementShellReissner4>();
            auto nodeid_vect = std::get<1>(el_it->second);
            std::array<std::shared_ptr<ChNodeFEAxyzrot>, 4> nodes;
            for (auto node_sel = 0; node_sel<4; ++node_sel)
            {
                // check if the node specified by the current element exists
                auto node = nodes_map.find(nodeid_vect[node_sel]);
                if (node != nodes_map.end())
                {
                    // check if the node specified by the current has not been inserted yet
                    auto node_found = inserted_nodes.find(nodeid_vect[node_sel]);
                    if (node_found == inserted_nodes.end())
                    {
                        nodes[node_sel] = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector<>(node->second[0], node->second[1], node->second[2]), QUNIT));
                        my_mesh->AddNode(nodes[node_sel]);
                        inserted_nodes.emplace_hint(inserted_nodes.end(), nodeid_vect[node_sel], nodes[node_sel]);
                    }
                    else
                    {
                        nodes[node_sel] = node_found->second;
                    }
                }
                else
                    throw ChException("Node not found\n");
            }
            // Add new element
            new_elem->SetNodes(nodes[3], nodes[2], nodes[1], nodes[0]);
            new_elem->AddLayer(shell_thickness, 0 * CH_C_DEG_TO_RAD, shell_material);
            //new_elem->SetAlphaDamp(0.0);
            my_mesh->AddElement(new_elem);
        }
    }

    GetLog() << "Added " << inserted_nodes.size() << " nodes over " << nodes_map.size() << ".\n";
    GetLog() << "Added " << my_mesh->GetElements().size() << " elements over " << elements_map.size() << ".\n";

    double external_threshold = 80.22e-3 - 2e-3;
    double internal_threshold = 25.5e-3 + 2e-3;
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> external_nodes;
    std::vector<std::shared_ptr<ChNodeFEAxyzrot>> internal_nodes;
    auto nodesmesh = my_mesh->GetNodes();
    for (auto node_sel = 0; node_sel<nodesmesh.size(); ++node_sel)
    {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(nodesmesh[node_sel]);
        if (node->GetPos().Length()>external_threshold)
        {
            external_nodes.push_back(node);
        }
        else if (node->GetPos().Length()<internal_threshold)
        {
            internal_nodes.push_back(node);
        }
    }

    for (auto it = internal_nodes.begin(); it!=internal_nodes.end(); ++it)
    {
        (*it)->SetFixed(true);
    }

    for (auto it = external_nodes.begin(); it != external_nodes.end(); ++it)
    {
        //(*it)->SetForce(ChVector<>(10,10,0));
        (*it)->SetForce(VECT_Z*100);
    }

    GetLog() << "External nodes: " << external_nodes.size() << "\n";
    GetLog() << "Internal nodes: " << internal_nodes.size() << "\n";



    //auto mvisualizeshellA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    //mvisualizeshellA->SetSmoothFaces(true);
    //mvisualizeshellA->SetWireframe(true);
    //my_mesh->AddAsset(mvisualizeshellA);

    
    auto mvisualizeshellB = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    //mvisualizeshellB->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizeshellB->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_STRESS_VONMISES);
    //mvisualizeshellB->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRESS);
    //mvisualizeshellB->SetWireframe(true);
    mvisualizeshellB->SetSymbolsThickness(0.01);
    my_mesh->AddAsset(mvisualizeshellB);
    

    //auto mvisualizeshellC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    //mvisualizeshellC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    // //mvisualizeshellC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_CSYS);
    //mvisualizeshellC->SetSymbolsThickness(0.05);
    //mvisualizeshellC->SetZbufferHide(false);
    //my_mesh->AddAsset(mvisualizeshellC);

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
    // Change solver to MKL
    auto mkl_solver = std::make_shared<ChSolverMKL<>>();
    mkl_solver->SetSparsityPatternLock(true);
    mkl_solver->ForceSparsityPatternUpdate(true);
    my_system.SetSolver(mkl_solver);

    /*
    my_system.SetSolverType(ChSolver::Type::MINRES); // <- NEEDED THIS or Matlab or MKL solver
    my_system.SetSolverWarmStarting(true); // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(200);
    my_system.SetMaxItersSolverStab(200);
    my_system.SetTolForce(1e-13);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);
    */

    // Change type of integrator:
    my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);
    // my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    // my_system.SetTimestepperType(ChTimestepper::NEWMARK);

    if (auto mint = std::dynamic_pointer_cast<ChImplicitIterativeTimestepper>(my_system.GetTimestepper())) {
        mint->SetMaxiters(5);
        mint->SetAbsTolerances(1e-12, 1e-12);
    }

    double timestep = 0.01;
    application.SetTimestep(timestep);
    my_system.Setup();
    my_system.Update();


    double mtime = 0;
    application.SetPaused(true);
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        // .. draw also a grid
        //ChIrrTools::drawGrid(application.GetVideoDriver(), 0.05, 0.05);

        if (!application.GetPaused())
        {
             //application.DoStep();
            // mtime = my_system.GetChTime();
            application.GetSystem()->DoStaticNonlinear(3);
            // application.GetSystem()->DoStaticLinear();
            mtime += timestep;
            GetLog() << "Update\n";
        }



        application.EndScene();

    }



    return 0;
}
