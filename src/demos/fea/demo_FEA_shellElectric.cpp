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
#include "chrono/physics/ChSystemNSC.h"


#include <set>
#include "chrono_fea/ChElementHexa_8.h"
#include "chrono_fea/ChElementShellReissner4.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChMeshFileLoader.h"
#include "chrono_fea/ChRotUtils.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_mkl/ChSolverMKL.h"
#include "chrono_postprocess/ChGnuPlot.h"
#include "utils/ChUtilsInputOutput.h"

// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;
using namespace chrono::postprocess;
using namespace irr;

// Output directory
std::string filename_sigma_t = "D:/SVN_MeltingLab/structural_EM/mesh/stressPriusCPSR.sigma_t.txt";
std::string filename_sigma_n = "D:/SVN_MeltingLab/structural_EM/mesh/stressPriusCPSR.sigma_n.txt";
//std::string filename_mesh = "D:/SVN_MeltingLab/structural_EM/mesh/prius_full_rotor_3D_5mm.INP";
std::string filename_mesh = "D:/SVN_MeltingLab/structural_EM/mesh/prius_3D_thickness_5mm_coarse.INP";
int test_num = 12345;
double omega = 5000*CH_C_2PI/60.0;

#define USE_MKL
#define USE_IRRLICHT
#define FULL_STRESS_OUTPUT
//#define EQUAL_ELEMENT_SPACING


class CSVwriter
{
private:
    std::ofstream myfile;
    std::ostringstream outbuffer;
public:

    CSVwriter(std::string filename)
    {
        myfile.open(filename, std::ios_base::out);
        if (!myfile.good())
            throw ChException("File with name: " + filename + "cannot be found");
    }

    template<typename T, typename... args_t>
    void AppendRow(const T& objects, const args_t&... objects_other)
    {
        outbuffer << objects << ", ";
        this->AppendRow(objects_other...);
    }

    template<typename T>
    void AppendRow(const T& objects)
    {
        outbuffer << objects;
        myfile << outbuffer.str() << std::endl;
        outbuffer.str("");
        outbuffer.clear();
    }


    ~CSVwriter()
    {
        myfile.close();
    }
};


void getSigmaGlob(const std::vector<ChVector<>>& sigma_glob_set, ChVector<>& sigma_glob, double angle)
{
    angle = fmod(angle, CH_C_2PI);
    if (angle < 0)
        angle += CH_C_2PI;

    assert(angle <= CH_C_2PI);
    assert(angle >= 0);

    double index = angle / CH_C_2PI * (sigma_glob_set.size()-1);
    auto index_int = static_cast<size_t>(floor(index));
    //int next_index = index_int > sigma_glob_set.size() ? index_int - sigma_glob_set.size() : index_int;
    int next_index = std::min(index_int + 1, sigma_glob_set.size()-1);

    sigma_glob = sigma_glob_set[index_int] + (index - index_int) * (sigma_glob_set[next_index] - sigma_glob_set[index_int]);
    
}

void getArcLength(const std::set<double>& angles, double& length_previous, double& length_next, double angle, double radius)
{
    angle = fmod(angle, CH_C_2PI);
    if (angle < 0)
        angle += CH_C_2PI;

    assert(angle <= CH_C_2PI);
    assert(angle >= 0);

    auto iter_prev = angles.lower_bound(angle);
    auto iter_next = iter_prev;
    double angle_previous = iter_prev != angles.begin() ? *--iter_prev : *(--angles.end()) - CH_C_2PI;
    double angle_next = (iter_next != angles.end()) && ++iter_next != angles.end() ? *iter_next : *(angles.begin())+CH_C_2PI;

    length_previous = 0.5*radius*(angle - angle_previous);
    length_next = 0.5*radius*(angle_next - angle);

    // OVERRIDE
    //length_next = 0.5*CH_C_2PI * 80.22e-3 / (560.0);
    //length_previous = length_next;

    assert(length_previous >= 0.0);
    assert(length_next >= 0.0);

}


int main(int argc, char* argv[]) {



    ChTimer<> tim;
    // Create a Chrono::Engine physical system
    ChSystemNSC my_system;

    auto truss = std::make_shared<ChBodyEasyBox>(0.1,0.1,0.1,7850);
    auto ass = std::make_shared<ChColorAsset>(1.0, 0.0, 0.0);
    truss->SetPos(ChVector<>(0.0, 0.0, -1.0));
    truss->AddAsset(ass);
    my_system.Add(truss);
    truss->SetBodyFixed(true);

#ifdef USE_IRRLICHT
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
                                   0.0, 0.5, 0, 512, video::SColorf((f32)0.8, (f32)0.8, (f32)1.0));

#endif

    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    auto my_mesh = std::make_shared<ChMesh>();

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);

    // my_system.Set_G_acc(VNULL); or
    my_mesh->SetAutomaticGravity(false);

    std::string element_tag = "C3D8";

    // Create a material
    double rho = 7850;
    double E = 200e9;
    double nu = 0.3;
    auto element_thickness = 5e-3;

    auto element_material = std::make_shared<ChContinuumElastic>(E, nu, rho);

    tim.start();
    std::map<std::shared_ptr<ChElementHexa_8>, unsigned int> inserted_elements_ptr_to_ID;
    std::map<std::shared_ptr<ChNodeFEAxyz>, unsigned int> inserted_nodes_ptr_to_ID;
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

    std::vector<ChVector<>> sigma_glob;
    sigma_glob.resize(sigma_n.size());
    auto delta_angle = CH_C_2PI / sigma_glob.size();
    for(auto sigma_sel = 0; sigma_sel<sigma_n.size(); ++sigma_sel)
    {

        sigma_glob[sigma_sel].x() = sigma_n[sigma_sel] * cos(delta_angle*sigma_sel) - sigma_t[sigma_sel] * sin(delta_angle*sigma_sel);
        sigma_glob[sigma_sel].y() = sigma_n[sigma_sel] * sin(delta_angle*sigma_sel) + sigma_t[sigma_sel] * cos(delta_angle*sigma_sel);
        sigma_glob[sigma_sel].z() = 0.0;
    }

    std::map<unsigned, std::tuple<std::string, std::vector<unsigned>>> elements_map;
    std::map<unsigned, std::vector<double>> nodes_map;
    std::map<std::string, std::vector<unsigned int>> nset_map;
    std::map<std::string, std::vector<unsigned int>> elset_map;
    std::map<unsigned int, std::shared_ptr<ChNodeFEAxyz>> inserted_nodes;

    try {
        ChMeshFileLoader::FromAbaqusFileMOD(filename_mesh, elements_map, nodes_map, nset_map, elset_map);
    }
    catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    // bool full_rotor = true;
    // int repetitions = 8;
    // std::list<std::shared_ptr<ChNodeFEAxyzrot>> nodes_to_check;

    // if (full_rotor) {
    //    for (auto slot_sel = 0; slot_sel < repetitions; ++slot_sel) {
    for (auto el_it = elements_map.begin(); el_it != elements_map.end(); ++el_it) {
        if (std::get<0>(el_it->second) == element_tag) {
            auto new_elem = std::make_shared<ChElementHexa_8>();
            auto nodeid_vect = std::get<1>(el_it->second);
            std::array<std::shared_ptr<ChNodeFEAxyz>, 8> nodes;
            for (auto node_sel = 0; node_sel < 8; ++node_sel) {
                // check if the node specified by the current element exists
                auto node = nodes_map.find(nodeid_vect[node_sel]);
                if (node != nodes_map.end()) {
                    // check if the node specified by the current has not been inserted yet
                    auto node_found = inserted_nodes.find(nodeid_vect[node_sel]);
                    if (node_found == inserted_nodes.end()) {
                        // nodes[node_sel] = std::make_shared<ChNodeFEAxyz>(ChFrame<>(ChVector<>(node->second[0],
                        // node->second[1], node->second[2]), QUNIT));  nodes[node_sel] =
                        // std::make_shared<ChNodeFEAxyz>(ChVector<>(node->second[0], node->second[1], node->second[2]),
                        // VECT_Z);
                        nodes[node_sel] =
                            std::make_shared<ChNodeFEAxyz>(ChVector<>(node->second[0], node->second[1], node->second[2]));
                        my_mesh->AddNode(nodes[node_sel]);
                        inserted_nodes.emplace_hint(inserted_nodes.end(), nodeid_vect[node_sel], nodes[node_sel]);
                        inserted_nodes_ptr_to_ID.emplace_hint(inserted_nodes_ptr_to_ID.end(), nodes[node_sel], nodeid_vect[node_sel]);
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
            
            my_mesh->AddElement(new_elem);
            inserted_elements_ptr_to_ID.emplace_hint(inserted_elements_ptr_to_ID.end(), new_elem, el_it->first);
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
    double external_threshold = rotor_external_radius - 1e-4;
    double internal_threshold = rotor_internal_radius + 1e-4;
    auto nodesmesh = my_mesh->GetNodes();
    std::vector<std::shared_ptr<ChNodeFEAxyz>> external_nodes;
    std::vector<std::shared_ptr<ChNodeFEAxyz>> internal_nodes;
    std::set<double> angles;

    for (auto node_sel = 0; node_sel < nodesmesh.size(); ++node_sel) {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(nodesmesh[node_sel]);
        auto dist_from_center = sqrt(node->GetPos().x()*node->GetPos().x() + node->GetPos().y()*node->GetPos().y());
        if (dist_from_center > external_threshold) {
            external_nodes.push_back(node);
            double angle = atan2(node->GetPos().y(), node->GetPos().x());
            if (angle < 0)
                angle += CH_C_2PI;
            angles.insert(angle);
        }
        else if (dist_from_center < internal_threshold) {
            internal_nodes.push_back(node);
        }
    }

    // fix nodes
    for (auto it = internal_nodes.begin(); it != internal_nodes.end(); ++it) {
        (*it)->SetFixed(true);
    }

    CSVwriter forces("forces.txt");
    for (auto it = external_nodes.begin(); it != external_nodes.end(); ++it) {
        double angle = atan2((*it)->GetPos().y(), (*it)->GetPos().x());

        // assure angle between [0, 2*pi)
        angle = fmod(angle, CH_C_2PI);
        if (angle < 0)
            angle += CH_C_2PI;

        assert(angle <= CH_C_2PI);
        assert(angle >= 0);

#ifdef EQUAL_ELEMENT_SPACING
        double index = angle / CH_C_2PI * (sigma_n.size() - 1);
        auto index_int = static_cast<size_t>(floor(index));
        //int next_index = index_int > sigma_glob_set.size() ? index_int - sigma_glob_set.size() : index_int;
        int next_index = std::min(index_int + 1, sigma_n.size() - 1);

        ChVector<> sigma_loc;
        sigma_loc[0] = sigma_n[index_int];
        sigma_loc[1] = sigma_t[index_int];
        sigma_loc[0] += (index - index_int) * (sigma_n[next_index] - sigma_n[index_int]);
        sigma_loc[1] += (index - index_int) * (sigma_t[next_index] - sigma_t[index_int]);
        sigma_loc[2] = 0.0;

        //ChMatrix33<double> rot_mat;
        //rot_mat(0, 0) = cos(angle);
        //rot_mat(1, 0) = sin(angle);
        //rot_mat(0, 1) = -sin(angle);
        //rot_mat(1, 1) = cos(angle);
        //rot_mat(2, 2) = 1.0;
        //ChVector<> forces_glob = rot_mat * sigma_loc;


        ChVector<> forces_glob;
        forces_glob[0] = sigma_loc[0] * cos(angle) - sigma_loc[1] * sin(angle);
        forces_glob[1] = sigma_loc[0] * sin(angle) + sigma_loc[1] * cos(angle);
        forces_glob[2] = 0;
        forces_glob.Scale(CH_C_2PI * rotor_external_radius * element_thickness / (2.0*external_nodes.size()));
        forces.AppendRow(angle, sigma_loc[0], sigma_loc[1], sigma_loc[2], forces_glob[0], forces_glob[1], forces_glob[2]);


#else
        auto iter_prev = angles.lower_bound(angle);
        auto iter_next = iter_prev;
        double halfangle_previous = 0.5*((iter_prev != angles.begin() ? *--iter_prev : *(--angles.end()))+angle);
        double halfangle_next = 0.5*(((iter_next != angles.end()) && ++iter_next != angles.end() ? *iter_next : *(angles.begin()))+angle);
        ChVector<> sigma_glob_previous, sigma_glob_next, sigma_glob_center;

        getSigmaGlob(sigma_glob, sigma_glob_previous, halfangle_previous);
        getSigmaGlob(sigma_glob, sigma_glob_next, halfangle_next);
        getSigmaGlob(sigma_glob, sigma_glob_center, angle);
        double archLength_previous, archLength_next;
        getArcLength(angles, archLength_previous, archLength_next, angle, rotor_external_radius);

        ChVector<> forces_glob = 0.5*(archLength_previous*element_thickness*sigma_glob_previous + archLength_next * element_thickness*sigma_glob_next);
        ChVector<> forces_glob2 = 0.5*(archLength_previous*element_thickness*sigma_glob_center + archLength_next * element_thickness*sigma_glob_center);

        forces.AppendRow(angle, 0.0,0.0,0.0, forces_glob[0], forces_glob[1], forces_glob[2]);

#endif



        //GetLog() << "Angle " << angle * 180.0 / CH_C_PI << "\n Forces" << forces_glob << "\n";

        (*it)->SetForce(forces_glob);
    }

    GetLog() << "External nodes: " << external_nodes.size() << "\n";
    GetLog() << "Internal nodes: " << internal_nodes.size() << "\n";
    //    }
    //}
    tim.stop();
    GetLog() << "Load and set element: " << tim() << "\n";

#ifdef USE_IRRLICHT

     //auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
     //mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
     //mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
     //mvisualizemesh->SetShrinkElements(true, 0.85);
     //mvisualizemesh->SetSmoothFaces(true);
     //my_mesh->AddAsset(mvisualizemesh);

    auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(false);
    mvisualizemeshref->SetDrawInUndeformedReference(false);
    my_mesh->AddAsset(mvisualizemeshref);

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRESS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_STRESS_VONMISES);
    mvisualizemeshC->SetSymbolsThickness(element_thickness);
    mvisualizemeshC->SetColorscaleMinMax(0,250e6);
    my_mesh->AddAsset(mvisualizemeshC);

    application.AssetBindAll();
    application.AssetUpdateAll();
#endif

    tim.reset();
    tim.start();
    my_system.SetupInitial();

    CSVwriter b("centrifugal_forces.txt");
    for (auto el_it = my_mesh->GetElements().begin(); el_it != my_mesh->GetElements().end(); ++el_it)
    {
        auto el = std::dynamic_pointer_cast<ChElementHexa_8>(*el_it);
        ChVector<> mean_node;
        for (auto node_sel = 0; node_sel<8; ++node_sel)
        {
            mean_node += std::dynamic_pointer_cast<ChNodeFEAxyz>(el->GetNodeN(node_sel))->GetPos();
        }
        mean_node *= 1.0 / 8.0;
        auto dist_from_center = sqrt(mean_node.x()*mean_node.x() + mean_node.y()*mean_node.y());
        auto centrifugal_force = el->GetVolume()*element_material->Get_density()*omega*omega*dist_from_center;
        for (auto node_sel = 0; node_sel<8; ++node_sel)
        {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(el->GetNodeN(node_sel));
            auto old_force = node->GetForce();

            auto centrifugal_force_vector = ChVector<>(mean_node.x(), mean_node.y(), 0.0);
            centrifugal_force_vector.Normalize();
            centrifugal_force_vector *= centrifugal_force / 8.0;
            node->SetForce(old_force + centrifugal_force_vector);
            b.AppendRow(mean_node.x(), mean_node.y(), centrifugal_force_vector.x(), centrifugal_force_vector.y(), centrifugal_force_vector.z());
        }
    }

    

//
// THE SOFT-REAL-TIME CYCLE
//
// Change solver to MKL
#ifdef USE_MKL
    auto mkl_solver = std::make_shared<ChSolverMKL<>>();
    mkl_solver->SetSparsityPatternLock(true);
    mkl_solver->ForceSparsityPatternUpdate(true);
    my_system.SetSolver(mkl_solver);
#else
    my_system.SetSolverType(ChSolver::Type::MINRES);  // <- NEEDED THIS or Matlab or MKL solver
    my_system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetMaxItersSolverSpeed(200);
    my_system.SetMaxItersSolverStab(200);
    my_system.SetTolForce(1e-13);
    auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system.GetSolver());
    msolver->SetVerbose(false);
    msolver->SetDiagonalPreconditioning(true);
#endif

    my_system.Setup();
    my_system.Update();

    tim.stop();
    GetLog() << "Setup Initial time: " << tim() << "\n";

    tim.reset();
    tim.start();
    application.GetSystem()->DoStaticLinear();
    tim.stop();
    GetLog() << "Simulation time: " << tim() << "\n";


    // Export results
    {
        tim.reset();
        tim.start();
        std::ofstream myfile;
        std::ostringstream filename_export;
        filename_export << "stress_" << test_num << ".txt";
        myfile.open(filename_export.str());

        for (auto el_it = my_mesh->GetElements().begin(); el_it != my_mesh->GetElements().end(); ++el_it)
        {
            std::ostringstream line;
            auto el = std::dynamic_pointer_cast<ChElementHexa_8>(*el_it);
            auto stress = el->GetStress(0, 0, 0);

            line << inserted_elements_ptr_to_ID.at(el) << ", "
                 << stress.GetEquivalentVonMises() << ", " 
#ifdef FULL_STRESS_OUTPUT
                << stress.XX() << ", "
                << stress.YY() << ", "
                << stress.ZZ() << ", "
                << stress.XY() << ", "
                << stress.YZ() << ", "
                << stress.XZ() << ", "
#endif
                 << std::endl;
            myfile << line.str();

        }
        myfile.close();

        tim.stop();
        GetLog() << "Export time: " << tim() << "\n";
    }
    GetLog() << "Done\n";

    //GetLog() << "forces: " << external_nodes[0]->GetForce() << "\n";
    //GetLog() << "pos-pos0: " << external_nodes[0]->GetPos() - external_nodes[0]->GetX0() << "\n";
    //GetLog() << "stress: " << std::dynamic_pointer_cast<ChElementHexa_8>(my_mesh->GetElements()[50])->GetStress(0,0,0) << "\n";
    //GetLog() << "strain: " << std::dynamic_pointer_cast<ChElementHexa_8>(my_mesh->GetElements()[50])->GetStrain(0,0,0) << "\n";
    //GetLog() << "VonMises: " << std::dynamic_pointer_cast<ChElementHexa_8>(my_mesh->GetElements()[50])->GetStress(0, 0, 0).GetEquivalentVonMises() << "\n";

    //for (auto el_it = my_mesh->GetElements().begin(); el_it!= my_mesh->GetElements().end(); ++ el_it)
    //{
    //    std::cout << std::dynamic_pointer_cast<ChElementHexa_8>(*el_it)->GetStress(0, 0, 0).GetEquivalentVonMises() << std::endl;
    //}

#ifdef USE_IRRLICHT
    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        if (!application.GetPaused()) {

        }

        application.EndScene();
    }
#endif

    return 0;
}