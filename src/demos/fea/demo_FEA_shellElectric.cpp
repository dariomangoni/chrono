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

#include <set>
#include "chrono_fea/ChElementShellReissner4.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChMeshFileLoader.h"
#include "chrono_fea/ChRotUtils.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_mkl/ChSolverMKL.h"
#include "chrono_postprocess/ChGnuPlot.h"
#include "utils/ChUtilsInputOutput.h"
#include "chrono_fea/ChElementHexahedron.h"
#include "chrono_fea/ChElementHexa_8.h"

// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;
using namespace chrono::postprocess;
using namespace irr;

// Output directory
std::string filename_sigma_t = "D:/SVN_MeltingLab/structural_EM/MATLAB/stressPriusCPSR.sigma_t.txt";
std::string filename_sigma_n = "D:/SVN_MeltingLab/structural_EM/MATLAB/stressPriusCPSR.sigma_n.txt";
std::string filename_mesh = "D:/SVN_MeltingLab/structural_EM/mesh/prius_3D_thickness_5mm_coarse.INP";

#define USE_MKL

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
    application.AddTypicalCamera(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.2)),
        core::vector3dfCH(ChVector<>(0.0, 0.0, 0.0)));
    // application.SetContactsDrawMode(irr::ChIrrTools::CONTACT_DISTANCES);

    application.AddLightWithShadow(core::vector3dfCH(ChVector<>(0.0, 0.0, 0.5)), core::vector3df(0.0, 0.0, 0.0), 0.1,
        0.1, 0.5, 0, 512, video::SColorf((f32)0.8, (f32)0.8, (f32)1.0));

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // my_system.Set_G_acc(VNULL); or
    my_mesh->SetAutomaticGravity(false);


    std::vector<double> sigma_t;
    std::vector<double> sigma_n;
    {
        // acquire forces
        std::ifstream fin(filename_sigma_t);
        if (fin.good())
            GetLog() << "Parsing Abaqus INP file: " << filename_sigma_t << "\n";
        else
            throw ChException("ERROR opening Abaqus .inp file: " + std::string(filename_sigma_t) + "\n");

        std::string tmp;
        char delim = ',';  // Ddefine the delimiter to split by
        double val;

        while (std::getline(fin, tmp, delim)) {
            // Provide proper checks here for tmp like if empty
            // Also strip down symbols like !, ., ?, etc.
            // Finally push it.
            std::istringstream stoken(tmp);
            stoken >> val;
            sigma_t.push_back(val);
        }
    }

    {
        // acquire forces
        std::ifstream fin(filename_sigma_n);
        if (fin.good())
            GetLog() << "Parsing Abaqus INP file: " << filename_sigma_n << "\n";
        else
            throw ChException("ERROR opening Abaqus .inp file: " + std::string(filename_sigma_n) + "\n");

        std::string tmp;
        char delim = ',';  // Ddefine the delimiter to split by
        double val;

        while (std::getline(fin, tmp, delim)) {
            // Provide proper checks here for tmp like if empty
            // Also strip down symbols like !, ., ?, etc.
            // Finally push it.
            std::istringstream stoken(tmp);
            stoken >> val;
            sigma_n.push_back(val);
        }
    }

    typedef ChElementHexa_8 element_type;
    typedef ChNodeFEAxyz nodes_type;
    typedef ChContinuumElastic material_type;
    std::string element_tag = "C3D8";

    // Create a material
    double rho = 7850;
    double E = 200e9;
    double nu = 0.7;

    auto element_material = std::make_shared<material_type>(rho, E, nu);

    std::map<unsigned, std::tuple<std::string, std::vector<unsigned>>> elements_map;
    std::map<unsigned, std::vector<double>> nodes_map;
    std::map<std::string, std::vector<unsigned int>> nset_map;
    std::map<std::string, std::vector<unsigned int>> elset_map;
    std::map<unsigned int, std::shared_ptr<nodes_type>> inserted_nodes;

    try {
        ChMeshFileLoader::FromAbaqusFileMOD(filename_mesh, elements_map,
            nodes_map, nset_map, elset_map);
    }
    catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    auto shell_thickness = 1e-2;


    // bool full_rotor = true;
    // int repetitions = 8;
    // std::list<std::shared_ptr<ChNodeFEAxyzrot>> nodes_to_check;

    // if (full_rotor) {
    //    for (auto slot_sel = 0; slot_sel < repetitions; ++slot_sel) {
    for (auto el_it = elements_map.begin(); el_it != elements_map.end(); ++el_it) {
        if (std::get<0>(el_it->second) == element_tag) {
            auto new_elem = std::make_shared<element_type>();
            auto nodeid_vect = std::get<1>(el_it->second);
            std::array<std::shared_ptr<nodes_type>, 8> nodes;
            for (auto node_sel = 0; node_sel < 8; ++node_sel) {
                // check if the node specified by the current element exists
                auto node = nodes_map.find(nodeid_vect[node_sel]);
                if (node != nodes_map.end()) {
                    // check if the node specified by the current has not been inserted yet
                    auto node_found = inserted_nodes.find(nodeid_vect[node_sel]);
                    if (node_found == inserted_nodes.end()) {
                        //nodes[node_sel] = std::make_shared<nodes_type>(ChFrame<>(ChVector<>(node->second[0], node->second[1], node->second[2]), QUNIT));
                        //nodes[node_sel] = std::make_shared<nodes_type>(ChVector<>(node->second[0], node->second[1], node->second[2]), VECT_Z);
                        nodes[node_sel] = std::make_shared<nodes_type>(ChVector<>(node->second[0], node->second[1], node->second[2]));
                        my_mesh->AddNode(nodes[node_sel]);
                        inserted_nodes.emplace_hint(inserted_nodes.end(), nodeid_vect[node_sel], nodes[node_sel]);
                    }
                    else {
                        nodes[node_sel] = node_found->second;
                    }
                }
                else
                    throw ChException("Node not found\n");
            }
            // Add new element
            new_elem->SetNodes(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
            new_elem->SetMaterial(element_material);
            //new_elem->AddLayer(shell_thickness, 0 * CH_C_DEG_TO_RAD, element_material);
            // new_elem->SetAlphaDamp(0.0);
            my_mesh->AddElement(new_elem);
        }
    }

    GetLog() << "Added " << inserted_nodes.size() << " nodes over " << nodes_map.size() << ".\n";
    GetLog() << "Added " << my_mesh->GetElements().size() << " elements over " << elements_map.size() << ".\n";

    //// Clean duplicated nodes
    //// pick lateral nodes
    // double lateral_threshold = 2e-3;
    // auto nodesmesh = my_mesh->GetNodes();
    // for (auto node_sel = 0; node_sel < nodesmesh.size(); ++node_sel) {
    //    auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(nodesmesh[node_sel]);
    //    if (abs(atan2(node->GetPos().y(), node->GetPos().x()) - slot_sel * CH_C_2PI / repetitions) < lateral_threshold
    //    || abs(atan2(node->GetPos().y(), node->GetPos().x()) - slot_sel * CH_C_2PI / repetitions) < lateral_threshold)
    //    {
    //        nodes_to_check.push_back(node);
    //    }
    //}

    double rotor_external_radius = 80.22e-3;
    double rotor_internal_radius = 25.5e-3;
    double external_threshold = rotor_external_radius - 1e-3;
    double internal_threshold = rotor_internal_radius + 1e-3;
    auto nodesmesh = my_mesh->GetNodes();
    std::vector<std::shared_ptr<nodes_type>> external_nodes;
    std::vector<std::shared_ptr<nodes_type>> internal_nodes;
    for (auto node_sel = 0; node_sel < nodesmesh.size(); ++node_sel) {
        auto node = std::dynamic_pointer_cast<nodes_type>(nodesmesh[node_sel]);
        if (node->GetPos().Length() > external_threshold) {
            external_nodes.push_back(node);
        }
        else if (node->GetPos().Length() < internal_threshold) {
            internal_nodes.push_back(node);
        }
    }

    for (auto it = internal_nodes.begin(); it != internal_nodes.end(); ++it) {
        (*it)->SetFixed(true);
    }

    for (auto it = external_nodes.begin(); it != external_nodes.end(); ++it) {
        double angle = atan2((*it)->GetPos().x(), (*it)->GetPos().y());
        if (angle < 0)
            angle += CH_C_2PI;
        ChMatrix33<double> rot_mat;
        rot_mat(0, 0) = cos(angle);
        rot_mat(1, 0) = sin(angle);
        rot_mat(0, 1) = -sin(angle);
        rot_mat(1, 1) = cos(angle);
        rot_mat(2, 2) = 1;

        ChVector<> sigma_loc;
        double index = angle / CH_C_2PI * sigma_t.size();
        int index_int = floor(index);
        sigma_loc[0] = sigma_n[index_int];
        sigma_loc[1] = sigma_t[index_int];
        sigma_loc[0] += (index - index_int)*(sigma_n[index_int + 1] - sigma_n[index_int]);
        sigma_loc[1] += (index - index_int)*(sigma_t[index_int + 1] - sigma_t[index_int]);

        ChVector<> forces_glob = rot_mat * sigma_loc;
        forces_glob.Scale(CH_C_2PI*rotor_external_radius / external_nodes.size()*shell_thickness);

        //GetLog() << "Angle " << angle << "\n Forces" << forces_glob << "\n";


        (*it)->SetForce(forces_glob);
    }

    GetLog() << "External nodes: " << external_nodes.size() << "\n";
    GetLog() << "Internal nodes: " << internal_nodes.size() << "\n";
    //    }
    //}





    //auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    //mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    //mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
    //mvisualizemesh->SetShrinkElements(true, 0.85);
    //mvisualizemesh->SetSmoothFaces(true);
    //my_mesh->AddAsset(mvisualizemesh);

    /*auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(true);
    mvisualizemeshref->SetDrawInUndeformedReference(true);
    my_mesh->AddAsset(mvisualizemeshref);*/

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRESS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_STRESS_VONMISES);
    mvisualizemeshC->SetSymbolsThickness(shell_thickness);
    my_mesh->AddAsset(mvisualizemeshC);

    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!

    application.AssetUpdateAll();

    // Mark completion of system construction
    my_system.SetupInitial();

    //
    // THE SOFT-REAL-TIME CYCLE
    //
     //Change solver to MKL
#ifdef USE_MKL
    auto mkl_solver = std::make_shared<ChSolverMKL<>>();
    mkl_solver->SetSparsityPatternLock(true);
    mkl_solver->ForceSparsityPatternUpdate(true);
    my_system.SetSolver(mkl_solver);
#else
    my_system.SetSolverType(ChSolver::Type::MINRES); // <- NEEDED THIS or Matlab or MKL solver
    my_system.SetSolverWarmStarting(true); // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(200);
    my_system.SetMaxItersSolverStab(200);
    my_system.SetTolForce(1e-13);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);
#endif

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
    //application.SetPaused(true);
    application.GetSystem()->DoStaticLinear();
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        // .. draw also a grid
        // ChIrrTools::drawGrid(application.GetVideoDriver(), 0.05, 0.05);

        if (!application.GetPaused()) {
            // application.DoStep();
            // mtime = my_system.GetChTime();
            //application.GetSystem()->DoStaticNonlinear(3);
             //application.GetSystem()->DoStaticLinear();

            /*auto nodesmesh = my_mesh->GetNodes();
            for (auto node_sel = 0; node_sel < nodesmesh.size(); ++node_sel) {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzrot>(nodesmesh[node_sel]);
            node->SetForce(node->GetForce()*1.0); // Is this only done for visualization purposes???
            }*/

            mtime += timestep;
            GetLog() << "Update\n";
        }

        application.EndScene();
    }

    return 0;
}