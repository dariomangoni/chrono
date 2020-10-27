//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be 
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

//
//   Demo code about 
//
//     - FEA for 3D beams



// Include some headers used by this tutorial...

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverPMINRES.h"
#include "chrono/solver/ChIterativeSolverLS.h"
//#include "chrono_pardisomkl/ChSolverPardisoMKL.h"


#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChVisualizationFEAmesh.h"
#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChLinkDirFrame.h"

#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_modelreduction/ChEigenAnalysis.h"
#include <typeinfo>

#include <unsupported/Eigen/SparseExtra>

#include <set>


#include "chrono/physics/ChLoaderUV.h"
#include "chrono/physics/ChLoadContainer.h"

#include "chrono/fea/ChElementTetra_4.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include <chrono\physics\ChLinkMotorRotationSpeed.h>
#include <chrono_pardisomkl\ChSolverPardisoMKL.h>

// Remember to use the namespace 'chrono' because all classes 
// of Chrono::Engine belong to this namespace and its children...



using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

bool enable_eigenanalysis = true;
bool save_matrix = false;


int main(int argc, char* argv[])
{
    // Create a Chrono::Engine physical system
    ChSystemSMC my_system;
    ChTimer<> m_timer_do_step;

    // Create the Irrlicht visualization (open the Irrlicht device, 
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Engine modes", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(-0.1f, 0.2f, -0.2f));

    // DeadBody object
    auto my_deadbody = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_deadbody);
    my_deadbody->SetBodyFixed(true);
    my_deadbody->SetName("DeadBody");


    //    // FINITE ELEMENT MESH
    // Create a mesh, that is a container for groups
    // of FEA elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each solid element in the mesh,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->Set_E(200e9);  // rubber 0.01e9, steel 200e9
    mmaterial->Set_v(0.3);
    mmaterial->Set_RayleighDampingK(0.00400);
    mmaterial->Set_density(7800);

    // Load an ABAQUS .INP tetrahedron mesh file from disk, defining a tetrahedron mesh.
    // Note that not all features of INP files are supported. Also, quadratic tetrahedrons are promoted to linear.
    // This is much easier thanfrcreating all nodes and elements via C++ programming.
    // Ex. you can generate these .INP files using Abaqus or exporting from the SolidWorks simulation tool.


    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets_rod;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("modelreduction/Rod_engine.INP").c_str(), mmaterial, node_sets_rod);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets_crank;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("modelreduction/Crank_engine.INP").c_str(), mmaterial, node_sets_crank);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets_piston;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("modelreduction/Piston_engine.INP").c_str(), mmaterial,
                                         node_sets_piston, ChVector<>(0, 0.1, 0),
                                         ChMatrix33<>(90 * CH_C_DEG_TO_RAD, VECT_Y));
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }

    // creazione corpi rigidi per vincoli

    auto my_crank_axle = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_crank_axle);
    my_crank_axle->SetName("crank_axle");

    auto my_crank_rod = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_crank_rod);
    my_crank_rod->SetName("crank_rod");

    auto my_rod_foot = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_rod_foot);
    my_rod_foot->SetName("rod_foot");

    auto my_rod_head = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_rod_head);
    my_rod_head->SetName("rod_head");

    auto my_piston_rod = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_piston_rod);
    my_piston_rod->SetName("piston_rod");

    auto my_piston_axle = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_piston_axle);
    my_piston_axle->SetName("piston_axle");

    //Fix link between flexible and rigid bodies
    // crank
    auto nodeset_crank_axle1 = "BC_NSET_CRANK_AXLE1";

    for (auto i = 0; i < node_sets_crank.at(nodeset_crank_axle1).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_crank[nodeset_crank_axle1][i]),
                          my_crank_axle);
        my_system.Add(mlink);
    }

    auto nodeset_crank_axle2 = "BC_NSET_CRANK_AXLE2";

    for (auto i = 0; i < node_sets_crank.at(nodeset_crank_axle2).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_crank[nodeset_crank_axle2][i]),
                          my_crank_axle);
        my_system.Add(mlink);
    }

    auto nodeset_crank_rod = "BC_NSET_CRANK_ROD";

    for (auto i = 0; i < node_sets_crank.at(nodeset_crank_rod).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_crank[nodeset_crank_rod][i]),
                          my_crank_rod);
        my_system.Add(mlink);
    }
    // rod
    auto nodeset_rod_foot = "BC_NSET_ROD_FOOT";

    for (auto i = 0; i < node_sets_rod.at(nodeset_rod_foot).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_rod[nodeset_rod_foot][i]),
                          my_rod_foot);
        my_system.Add(mlink);
    }

    auto nodeset_rod_head = "BC_NSET_ROD_HEAD";

    for (auto i = 0; i < node_sets_rod.at(nodeset_rod_head).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_rod[nodeset_rod_head][i]), my_rod_head);
        my_system.Add(mlink);
    }
    //piston
    auto nodeset_piston_rod1 = "BC_NSET_PISTON_ROD1";

    for (auto i = 0; i < node_sets_piston.at(nodeset_piston_rod1).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_piston[nodeset_piston_rod1][i]),
                          my_piston_rod);
        my_system.Add(mlink);
    }

    auto nodeset_piston_rod2 = "BC_NSET_PISTON_ROD2";

    for (auto i = 0; i < node_sets_piston.at(nodeset_piston_rod2).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_piston[nodeset_piston_rod2][i]),
                          my_piston_rod);
        my_system.Add(mlink);
    }

    auto nodeset_piston_axle = "BC_NSET_PISTON_AXLE";

    for (auto i = 0; i < node_sets_piston.at(nodeset_piston_axle).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets_piston[nodeset_piston_axle][i]),
                          my_piston_axle);
        my_system.Add(mlink);
    }

    //// vincoli tra i corpi rigidi
    //// link between crank and reference frame
    //auto my_link_crank_axle = chrono_types::make_shared<ChLinkLockRevolute>();
    //my_link_crank_axle->SetName("RevJointCrankAxle");
    //my_link_crank_axle->Initialize(my_crank_axle, my_deadbody, ChCoordsys<>(ChVector<>(2, 0, 0)));
    //my_system.AddLink(my_link_crank_axle);


    auto my_link_AB = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    my_link_AB->Initialize(my_crank_axle, my_deadbody,
                           ChFrame<>(ChVector<>(0, 0, 0), 90 * CH_C_DEG_TO_RAD, ChVector<>(0, 1, 0)));
    my_link_AB->SetName("RotationalMotor");
    my_system.AddLink(my_link_AB);
    auto my_speed_function = chrono_types::make_shared<ChFunction_Const>(10*CH_C_PI);  // speed w=3.145 rad/sec
    my_link_AB->SetSpeedFunction(my_speed_function);

    //vincolo biella manovella
    auto my_link_BC = chrono_types::make_shared<ChLinkLockRevolute>();
    my_link_BC->SetName("RevJointCrankRod");
    my_link_BC->Initialize(my_crank_rod, my_rod_foot, ChCoordsys<>(ChVector<>(0, 0.035, 0), 90*CH_C_DEG_TO_RAD,VECT_Y));
    my_system.AddLink(my_link_BC);

    //vincolo biella pistone
    auto my_link_CD = chrono_types::make_shared<ChLinkLockRevolute>();
    my_link_CD->SetName("RevJointPistonRod");
    my_link_CD->Initialize(my_rod_head, my_piston_rod,
                           ChCoordsys<>(ChVector<>(0, 0.125, 0), 90 * CH_C_DEG_TO_RAD, VECT_Y));
    my_system.AddLink(my_link_CD);

    // vincolo pistone
    auto my_link_CA = chrono_types::make_shared<ChLinkLockCylindrical>();
    my_link_CA->SetName("TransJointRodGround");
    my_link_CA->Initialize(my_piston_axle, my_deadbody, ChCoordsys<>(ChVector<>(0, 0, 0), -90 * CH_C_DEG_TO_RAD, VECT_X));
    my_system.AddLink(my_link_CA);


    // We do not want gravity effect on FEA elements in this demo
    my_mesh->SetAutomaticGravity(false);

    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);


    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChTriangleMeshShape
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.
    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    // postprocessor that can handle a coloured ChTriangleMeshShape).
    // Do not forget AddAsset() at the end!


    /*
    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizebeamA->SetSmoothFaces(true);
    my_mesh->AddAsset(mvisualizebeamA);
    */

    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    //mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
    //mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizebeamA->SetDefaultMeshColor(ChColor(0,1,0,0));
    mvisualizebeamA->SetSmoothFaces(true);
    mvisualizebeamA->SetWireframe(false);
    my_mesh->AddAsset(mvisualizebeamA);

    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*my_mesh.get());
    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_CSYS);
    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizebeamC->SetSymbolsThickness(0.006);
    mvisualizebeamC->SetSymbolsScale(0.01);
    mvisualizebeamC->SetZbufferHide(false);
    my_mesh->AddAsset(mvisualizebeamC);


    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
    // in the system. These ChIrrNodeAsset assets are 'proxies' to the Irrlicht meshes.
    // If you need a finer control on which item really needs a visualization proxy in 
    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.

    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!

    application.AssetUpdateAll();

    //// Mark completion of system construction
    //my_system.SetupInitial();

    //auto solver = chrono_types::make_shared<ChSolverMINRES>();
    //my_system.SetSolver(solver);
    //solver->SetMaxIterations(600);
    //solver->SetTolerance(1e-8);
    //solver->EnableDiagonalPreconditioner(true);
    //solver->SetVerbose(true);

    // Configure PardisoMKL solver.
    // For this simple and relatively small problem, use of the sparsity pattern learner may introduce additional
    // overhead (if the sparsity pattern is not locked).
    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    mkl_solver->UseSparsityPatternLearner(true);
    mkl_solver->LockSparsityPattern(false);
    //mkl_solver->SetVerbose(false);
    my_system.SetSolver(mkl_solver);

    application.SetTimestep(0.001);

    std::cout << "finito" << std::endl;

    m_timer_do_step.start();
    application.DoStep();
    m_timer_do_step.stop();
    std::cout << "Time for do step:" << m_timer_do_step() << std::endl;

    // Save matrix
    if (save_matrix) {
    
    ChSparseMatrix matKaug;
    application.GetSystem()->KRMmatricesLoad(1.0, 0, 0);
    application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
    application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matKaug, nullptr, true);
    //std::cout << matKaug << std::endl;
    Eigen::saveMarket(matKaug, "matKaug.mat");

    ChSparseMatrix matMaug;
    application.GetSystem()->KRMmatricesLoad(0, 0, 1.0);
    application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
    application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matMaug, nullptr, false);
    //std::cout << matMaug << std::endl;
    Eigen::saveMarket(matMaug, "matMaug.mat");

    ChSparseMatrix matM_expl;
    application.GetSystem()->GetMassMatrix(&matM_expl);
    matM_expl.makeCompressed();
    //std::cout << matM_expl << std::endl;
    Eigen::saveMarket(matM_expl, "matM_expl.mat");

    ChSparseMatrix matK_expl;
    application.GetSystem()->GetStiffnessMatrix(&matK_expl);
    matK_expl.makeCompressed();
    //std::cout << matK_expl << std::endl;
    Eigen::saveMarket(matK_expl, "matK_expl.mat");

    ChSparseMatrix matCq_expl;
    application.GetSystem()->GetConstraintJacobianMatrix(&matCq_expl);
    matCq_expl.makeCompressed();
    //std::cout << matCq_expl << std::endl;
    Eigen::saveMarket(matCq_expl, "matCq_expl.mat");
    }


    if (enable_eigenanalysis) {
        GetLog() << "\n\n===========EIGENPROBLEM======== \n";

        std::cout << "Starting EigenSolve" << std::endl;

        ChEigenAnalysis eig_analysis(application);
        eig_analysis.SetVerbose(true);
        eig_analysis.EigenAnalysis(10, 1e4, false);

        std::cout << "Ended EigenSolve" << std::endl;

        ChVectorDynamic<> residualsNorm;
        eig_analysis.GetResidualsNorm(residualsNorm);
        std::cout << "Residual norm rev: " << residualsNorm.reverse() << std::endl;

        ChVectorDynamic<> frequencies;
        eig_analysis.GetFrequencies(frequencies);
        std::cout << "frequencies rev: " << frequencies.reverse() << std::endl;

        std::cout << "Time for assembly:" << eig_analysis.GetTime_Assembly() << std::endl;
        std::cout << "Time for the eigensolver:" << eig_analysis.GetTime_Eigensolve() << std::endl;


        Eigen::saveMarket(eig_analysis.GetEigenValues(), "eig_val.mat");
        Eigen::saveMarket(eig_analysis.GetEigenVectors(), "eig_vect.mat");


        while (application.GetDevice()->run()) {
            application.BeginScene();
            application.DrawAll();
            eig_analysis.UpdateEigenMode();
            application.EndScene();
        }
        
    }
    else {
        while (application.GetDevice()->run()) {
            application.BeginScene();
            application.DrawAll();
            application.DoStep();
            application.EndScene();
        }
    }

}


