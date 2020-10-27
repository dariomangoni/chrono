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


    // Create the Irrlicht visualization (open the Irrlicht device, 
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"Rod modes", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(-0.1f, 0.2f, -0.2f));


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


    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets;

    try {
        ChMeshFileLoader::FromAbaqusFile(my_mesh, GetChronoDataFile("modelreduction/Rod.INP").c_str() , mmaterial, node_sets);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }
    //


    // DeadBody object
    auto my_deadbody = chrono_types::make_shared<ChBody>();
    my_system.AddBody(my_deadbody);
    my_deadbody->SetBodyFixed(true);
    my_deadbody->SetName("DeadBody");

    //Fix link between beam and deadboy

    auto nodeset_sel1 = "BC_NSET1";
    for (auto i = 0; i < node_sets.at(nodeset_sel1).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets[nodeset_sel1][i]), my_deadbody);
        my_system.Add(mlink);
    }

    auto nodeset_sel2 = "BC_NSET2";
    for (auto i = 0; i < node_sets.at(nodeset_sel2).size(); ++i) {
        auto mlink = chrono_types::make_shared<ChLinkPointFrame>();
        mlink->Initialize(std::dynamic_pointer_cast<ChNodeFEAxyz>(node_sets[nodeset_sel2][i]), my_deadbody);
        my_system.Add(mlink);
    }

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

    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    my_system.SetSolver(solver);
    solver->SetMaxIterations(600);
    solver->SetTolerance(1e-13);
    solver->EnableDiagonalPreconditioner(true);
    solver->SetVerbose(true);

    // Configure PardisoMKL solver.
    // For this simple and relatively small problem, use of the sparsity pattern learner may introduce additional
    // overhead (if the sparsity pattern is not locked).
    //auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    //mkl_solver->UseSparsityPatternLearner(false);
    //mkl_solver->LockSparsityPattern(false);
    //mkl_solver->SetVerbose(false);
    //my_system.SetSolver(mkl_solver);

    application.SetTimestep(0.001);

    application.DoStep();

    // Save matrix
    if (save_matrix) {
        ChSparseMatrix matKaug;
        application.GetSystem()->KRMmatricesLoad(1.0, 0, 0);
        application.GetSystem()->GetSystemDescriptor()->SetMassFactor(0.0);
        application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matKaug, nullptr, true);
        // std::cout << matKaug << std::endl;
        Eigen::saveMarket(matKaug, "matKaug.mat");

        ChSparseMatrix matMaug;
        application.GetSystem()->KRMmatricesLoad(0, 0, 1.0);
        application.GetSystem()->GetSystemDescriptor()->SetMassFactor(1.0);
        application.GetSystem()->GetSystemDescriptor()->ConvertToMatrixForm(&matMaug, nullptr, false);
        // std::cout << matMaug << std::endl;
        Eigen::saveMarket(matMaug, "matMaug.mat");

        ChSparseMatrix matM_expl;
        application.GetSystem()->GetMassMatrix(&matM_expl);
        matM_expl.makeCompressed();
        // std::cout << matM_expl << std::endl;
        Eigen::saveMarket(matM_expl, "matM_expl.mat");

        ChSparseMatrix matK_expl;
        application.GetSystem()->GetStiffnessMatrix(&matK_expl);
        matK_expl.makeCompressed();
        // std::cout << matK_expl << std::endl;
        Eigen::saveMarket(matK_expl, "matK_expl.mat");

        ChSparseMatrix matCq_expl;
        application.GetSystem()->GetConstraintJacobianMatrix(&matCq_expl);
        matCq_expl.makeCompressed();
        // std::cout << matCq_expl << std::endl;
        Eigen::saveMarket(matCq_expl, "matCq_expl.mat");
    }

    if (enable_eigenanalysis) {
        GetLog() << "\n\n===========EIGENPROBLEM======== \n";

        ChEigenAnalysis eig_analysis(application);
        eig_analysis.SetVerbose(true);
        eig_analysis.EigenAnalysis(10, 10e7, false);

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

        // Strict realtime timer
        std::function<void()> simadvance_fun = std::bind(&ChEigenAnalysis::UpdateEigenMode, &eig_analysis, 1);
        // std::function<void()> render_fun = std::bind(&ChIrrApp::DrawAll, application);
        std::function<void(double)> render_fun = [&application](double dummy) { application.DrawAll(); };
        std::function<void()> pre_sim_fun = std::bind(&ChIrrApp::BeginScene, &application, true, true, irr::video::SColor(255, 0, 0, 0));

        // std::function<void(double)> render_fun = [&application, &eig_analysis](double ratio) {
        // eig_analysis.UpdateEigenMode(ratio); application.DrawAll(); };
        ChRealtimeDualStepTimer drealtime = {application, simadvance_fun, render_fun};
        drealtime.SetPreAdvanceSimulationFunction(pre_sim_fun);

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


