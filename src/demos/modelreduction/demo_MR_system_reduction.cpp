#include "chrono_modelreduction/ChModelReduction.h"
#include <Eigen/Eigenvalues> 
////
//#include <Eigen/Core>
//#include <Eigen/SparseCore>
//#include <Eigen/Eigenvalues>
#include <SymGEigsSolver.h>
//#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseCholesky.h>
//#include <iostream>
#include <MatOp/SparseSymMatProd.h>
#include "core/ChCSR3Matrix.h"
using namespace Spectra;
////
//


//
//
//int test_eigen()
//{
//    const int n = 10;
//    chrono::ChCSR3Matrix mymatA(n, n);
//
//    mymatA.SetElement(0, 0, 1);
//    mymatA.SetElement(1, 1, 2);
//    mymatA.SetElement(2, 2, 3);
//    mymatA.SetElement(3, 3, 4);
//    mymatA.SetElement(4, 4, 5);
//    mymatA.SetElement(5, 5, 6);
//    mymatA.SetElement(6, 6, 7);
//    mymatA.SetElement(7, 7, 8);
//    mymatA.SetElement(8, 8, 9);
//    mymatA.SetElement(9, 9, 10);
//    mymatA.Compress();
//
//    Eigen::SparseMatrix<double, Eigen::RowMajor> mymatA_wrap(n, n);
//    mymatA_wrap.LoadChCSR3Matrix(mymatA);
//
//    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
//    chrono::ChCSR3Matrix mymatB(n, n);
//
//    mymatB.SetElement(0, 0, 1);
//    mymatB.SetElement(1, 1, 2);
//    mymatB.SetElement(2, 2, 3);
//    mymatB.SetElement(3, 3, 4);
//    mymatB.SetElement(4, 4, 5);
//    mymatB.SetElement(5, 5, 6);
//    mymatB.SetElement(6, 6, 7);
//    mymatB.SetElement(7, 7, 8);
//    mymatB.SetElement(8, 8, 9);
//    mymatB.SetElement(9, 9, 10);
//    mymatB.Compress();
//
//
//    Eigen::SparseMatrix<double, Eigen::RowMajor> mymatB_wrap(n, n);
//    mymatB_wrap.LoadChCSR3Matrix(mymatB);
//
//    // Construct matrix operation object using the wrapper classes
//    SparseSymMatProd<double> op(mymatA_wrap);
//    SparseCholesky<double>  Bop(mymatB_wrap);
//
//    //// Verify results using the generalized eigen solver in Eigen
//    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA_wrap, mymatB_wrap);
//    //Eigen::GeneralizedEigenSolver<Eigen::SparseMatrix<double>> es(mymatA_wrap, mymatB_wrap, true);
//    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
//    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
//    return 0;
//}
//
//
//
//int model_reduction_generalized_useroperation2()
//{
//    // We are going to solve the generalized eigenvalue problem A * x = lambda * B * x
//    const int n = 10;
//
//    // Define the A matrix
//    Eigen::SparseMatrix<double> A(n, n);
//    for (int i = 0; i < n; i++)
//    {
//        A.insert(i, i) = i + 1;
//    }
//
//    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
//    Eigen::SparseMatrix<double> B(n, n);
//    B.reserve(Eigen::VectorXi::Constant(n, 3));
//    for (int i = 0; i < n; i++)
//    {
//        B.insert(i, i) = 2.0;
//        if (i > 0)
//            B.insert(i - 1, i) = 1.0;
//        if (i < n - 1)
//            B.insert(i + 1, i) = 1.0;
//    }
//
//    // Construct matrix operation object using the wrapper classes
//    SparseSymMatProd<double> op(A);
//    SparseCholesky<double>  Bop(B);
//
//    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//    SymGEigsSolver<double, LARGEST_ALGE, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);
//
//    // Initialize and compute
//    geigs.init();
//    int nconv = geigs.compute();
//
//    // Retrieve results
//    Eigen::VectorXd evalues;
//    Eigen::MatrixXd evecs;
//    if (geigs.info() == SUCCESSFUL)
//    {
//        evalues = geigs.eigenvalues();
//        evecs = geigs.eigenvectors();
//    }
//    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
//    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;
//
//    return 0;
//}
//
//
//
//
//
//int main()
//{
//    model_reduction_generalized_useroperation();
//    return 0;
//}
//
//
//
//
//
//
//////
////// PROJECT CHRONO - http://projectchrono.org
//////
////// Copyright (c) 2013 Project Chrono
////// All rights reserved.
//////
////// Use of this source code is governed by a BSD-style license that can be 
////// found in the LICENSE file at the top level of the distribution
////// and at http://projectchrono.org/license-chrono.txt.
//////
////
//////
//////   Demo code about 
//////
//////     - FEA for 3D beams
////
////
////
////// Include some headers used by this tutorial...
////
////#include "chrono/physics/ChSystem.h"
////#include "chrono/physics/ChLinkMate.h"
////#include "chrono/physics/ChBodyEasy.h"
////#include "chrono/timestepper/ChTimestepper.h"
////#include "chrono/solver/ChSolverPMINRES.h"
////#include "chrono/solver/ChSolverMINRES.h"
////
////#include "chrono_fea/ChElementBeamEuler.h"
////#include "chrono_fea/ChBuilderBeam.h"
////#include "chrono_fea/ChMesh.h"
////#include "chrono_fea/ChVisualizationFEAmesh.h"
////#include "chrono_fea/ChLinkPointFrame.h"
////#include "chrono_fea/ChLinkDirFrame.h"
////
////#include "chrono_irrlicht/ChIrrApp.h"
////#include <chrono_modelreduction/ChModelReduction.h>
////
////// Remember to use the namespace 'chrono' because all classes 
////// of Chrono::Engine belong to this namespace and its children...
////
////using namespace chrono;
////using namespace chrono::fea;
////using namespace chrono::irrlicht;
////
////using namespace irr;
////
////
////
////int main(int argc, char* argv[])
////{
////
////    //model_reduction_usermultiplication();
////    //model_reduction_test2();
////    model_reduction_generalized2();
////
////    // Create a Chrono::Engine physical system
////    ChSystem my_system;
////
////
////    // Create the Irrlicht visualization (open the Irrlicht device, 
////    // bind a simple user interface, etc. etc.)
////    ChIrrApp application(&my_system, L"Beams (SPACE for dynamics, F10 / F11 statics)", core::dimension2d<u32>(800, 600), false, true);
////
////    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
////    application.AddTypicalLogo();
////    application.AddTypicalSky();
////    application.AddTypicalLights();
////    application.AddTypicalCamera(core::vector3df(-0.1f, 0.2f, -0.2f));
////
////
////
////    // Create a mesh, that is a container for groups
////    // of elements and their referenced nodes.
////    auto my_mesh = std::make_shared<ChMesh>();
////
////
////    // Create a section, i.e. thickness and material properties
////    // for beams. This will be shared among some beams.
////
////    auto msection = std::make_shared<ChBeamSectionAdvanced>();
////
////    double beam_wy = 0.012;
////    double beam_wz = 0.025;
////    msection->SetAsRectangularSection(beam_wy, beam_wz);
////    msection->SetYoungModulus(0.01e9);
////    msection->SetGshearModulus(0.01e9 * 0.3);
////    msection->SetBeamRaleyghDamping(0.000);
////    //msection->SetCentroid(0,0.02); 
////    //msection->SetShearCenter(0,0.1); 
////    //msection->SetSectionRotation(45*CH_C_RAD_TO_DEG);
////
////
////    //
////    // Add some EULER-BERNOULLI BEAMS:
////    //
////
////    double beam_L = 0.1;
////
////
////    //auto hnode1 = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector<>(0, 0, 0)));
////    //auto hnode2 = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector<>(beam_L, 0, 0)));
////    //auto hnode3 = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(ChVector<>(beam_L * 2, 0, 0)));
////
////    //my_mesh->AddNode(hnode1);
////    //my_mesh->AddNode(hnode2);
////    //my_mesh->AddNode(hnode3);
////
////    //auto belement1 = std::make_shared<ChElementBeamEuler>();
////
////    //belement1->SetNodes(hnode1, hnode2);
////    //belement1->SetSection(msection);
////
////    //my_mesh->AddElement(belement1);
////
////
////    //auto belement2 = std::make_shared<ChElementBeamEuler>();
////
////    //belement2->SetNodes(hnode2, hnode3);
////    //belement2->SetSection(msection);
////
////    //my_mesh->AddElement(belement2);
////
////
////    //// Apply a force or a torque to a node:
////    //hnode2->SetForce(ChVector<>(4, 2, 0));
////    ////hnode3->SetTorque( ChVector<>(0, -0.04, 0));
////
////
////    //// Fix a node to ground:
////    ////hnode1->SetFixed(true);
////    //auto mtruss = std::make_shared<ChBody>();
////    //mtruss->SetBodyFixed(true);
////    //my_system.Add(mtruss);
////
////    //auto constr_bc = std::make_shared<ChLinkMateGeneric>();
////    //constr_bc->Initialize(hnode3,
////    //                      mtruss,
////    //                      false,
////    //                      hnode3->Frame(),
////    //                      hnode3->Frame()
////    //);
////    //my_system.Add(constr_bc);
////    //constr_bc->SetConstrainedCoords(true, true, true,	  // x, y, z
////    //                                true, true, true);   // Rx, Ry, Rz
////
////    //auto constr_d = std::make_shared<ChLinkMateGeneric>();
////    //constr_d->Initialize(hnode1,
////    //                     mtruss,
////    //                     false,
////    //                     hnode1->Frame(),
////    //                     hnode1->Frame()
////    //);
////    //my_system.Add(constr_d);
////    //constr_d->SetConstrainedCoords(false, true, true,	  // x, y, z
////                                   //false, false, false);   // Rx, Ry, Rz
////
////
////    //
////    // Add some EULER-BERNOULLI BEAMS (the fast way!)
////    //
////
////    // Shortcut!
////    // This ChBuilderBeam helper object is very useful because it will 
////    // subdivide 'beams' into sequences of finite elements of beam type, ex.
////    // one 'beam' could be made of 5 FEM elements of ChElementBeamEuler class.
////    // If new nodes are needed, it will create them for you.
////    ChBuilderBeam builder;
////
////    // Now, simply use BuildBeam to create a beam from a point to another: 
////    builder.BuildBeam(my_mesh,		// the mesh where to put the created nodes and elements 
////                      msection,		// the ChBeamSectionAdvanced to use for the ChElementBeamEuler elements
////                      3,				// the number of ChElementBeamEuler to create
////                      ChVector<>(0, 0, -0.1),		// the 'A' point in space (beginning of beam)
////                      ChVector<>(0.2, 0, -0.1),	// the 'B' point in space (end of beam)
////                      ChVector<>(0, 1, 0));			// the 'Y' up direction of the section for the beam
////
////                                                    // After having used BuildBeam(), you can retrieve the nodes used for the beam,
////                                                    // For example say you want to fix the A end and apply a force to the B end:
////    builder.GetLastBeamNodes().back()->SetFixed(true);
////    builder.GetLastBeamNodes().front()->SetForce(ChVector<>(0, -1, 0));
////
////    //// Again, use BuildBeam for creating another beam, this time
////    //// it uses one node (the last node created by the last beam) and one point:
////    //builder.BuildBeam(my_mesh,
////    //                  msection,
////    //                  5,
////    //                  builder.GetLastBeamNodes().front(), // the 'A' node in space (beginning of beam)
////    //                  ChVector<>(0.2, 0.1, -0.1),	// the 'B' point in space (end of beam)
////    //                  ChVector<>(0, 1, 0));			// the 'Y' up direction of the section for the beam
////
////
////
////
////
////    //
////    // Final touches..
////    // 
////
////
////    // We do not want gravity effect on FEA elements in this demo
////    my_mesh->SetAutomaticGravity(false);
////
////    // Remember to add the mesh to the system!
////    my_system.Add(my_mesh);
////
////
////
////
////    // ==Asset== attach a visualization of the FEM mesh.
////    // This will automatically update a triangle mesh (a ChTriangleMeshShape
////    // asset that is internally managed) by setting  proper
////    // coordinates and vertex colours as in the FEM elements.
////    // Such triangle mesh can be rendered by Irrlicht or POVray or whatever
////    // postprocessor that can handle a coloured ChTriangleMeshShape).
////    // Do not forget AddAsset() at the end!
////
////
////    /*
////    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
////    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
////    mvisualizebeamA->SetSmoothFaces(true);
////    my_mesh->AddAsset(mvisualizebeamA);
////    */
////
////    auto mvisualizebeamA = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
////    mvisualizebeamA->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_ELEM_BEAM_MZ);
////    mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
////    mvisualizebeamA->SetSmoothFaces(true);
////    mvisualizebeamA->SetWireframe(false);
////    my_mesh->AddAsset(mvisualizebeamA);
////
////    auto mvisualizebeamC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
////    mvisualizebeamC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_CSYS);
////    mvisualizebeamC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
////    mvisualizebeamC->SetSymbolsThickness(0.006);
////    mvisualizebeamC->SetSymbolsScale(0.01);
////    mvisualizebeamC->SetZbufferHide(false);
////    my_mesh->AddAsset(mvisualizebeamC);
////
////
////    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
////    // in the system. These ChIrrNodeAsset assets are 'proxies' to the Irrlicht meshes.
////    // If you need a finer control on which item really needs a visualization proxy in 
////    // Irrlicht, just use application.AssetBind(myitem); on a per-item basis.
////
////    application.AssetBindAll();
////
////    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
////    // that you added to the bodies into 3D shapes, they can be visualized by Irrlicht!
////
////    application.AssetUpdateAll();
////
////    // Mark completion of system construction
////    my_system.SetupInitial();
////
////
////    // 
////    // THE SOFT-REAL-TIME CYCLE
////    //
////    my_system.SetSolverType(ChSystem::SOLVER_MINRES);
////    my_system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
////    my_system.SetMaxItersSolverSpeed(460);
////    my_system.SetMaxItersSolverStab(460);
////    my_system.SetTolForce(1e-13);
////    ChSolverMINRES* msolver = (ChSolverMINRES*)my_system.GetSolverSpeed();
////    msolver->SetVerbose(false);
////    msolver->SetDiagonalPreconditioning(true);
////
////    /*
////    // TEST: The Matlab external linear solver, for max precision in benchmarks
////    ChMatlabEngine matlab_engine;
////    ChMatlabSolver* matlab_solver_stab  = new ChMatlabSolver(matlab_engine);
////    ChMatlabSolver* matlab_solver_speed = new ChMatlabSolver(matlab_engine);
////    my_system.ChangeSolverStab (matlab_solver_stab);
////    my_system.ChangeSolverSpeed(matlab_solver_speed);
////    application.GetSystem()->Update();
////    application.SetPaused(true);
////    */
////
////
////    // Change type of integrator: 
////    my_system.SetIntegrationType(chrono::ChSystem::INT_HHT);
////
////    // if later you want to change integrator settings:
////    if (auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper())) {
////        mystepper->SetAlpha(-0.2);
////        mystepper->SetMaxiters(6);
////        mystepper->SetAbsTolerances(1e-12);
////        mystepper->SetVerbose(true);
////        mystepper->SetStepControl(false);
////    }
////
////    my_system.SetIntegrationType(chrono::ChSystem::INT_EULER_IMPLICIT_LINEARIZED);
////
////    application.SetTimestep(0.001);
////
////
////    GetLog() << "\n\n\n===========STATICS======== \n\n\n";
////
////
////
////    //	application.GetSystem()->DoStaticLinear();
////
////
////    GetLog() << "BEAM RESULTS (LINEAR STATIC ANALYSIS) \n\n";
////
////
////    //ChVector<> F, M;
////    //ChMatrixDynamic<> displ;
////
////    //belement1->GetStateBlock(displ);
////    //GetLog() << displ;
////    //for (double eta = -1; eta <= 1; eta += 0.4)
////    //{
////    //    belement1->EvaluateSectionForceTorque(eta, displ, F, M);
////    //    GetLog() << "  b1_at " << eta << " Mx=" << M.x << " My=" << M.y << " Mz=" << M.z << " Tx=" << F.x << " Ty=" << F.y << " Tz=" << F.z << "\n";
////    //}
////    //GetLog() << "\n";
////    //belement2->GetStateBlock(displ);
////    //for (double eta = -1; eta <= 1; eta += 0.4)
////    //{
////    //    belement2->EvaluateSectionForceTorque(eta, displ, F, M);
////    //    GetLog() << "  b2_at " << eta << " Mx=" << M.x << " My=" << M.y << " Mz=" << M.z << " Tx=" << F.x << " Ty=" << F.y << " Tz=" << F.z << "\n";
////    //}
////
////    //GetLog() << "Node 3 coordinate x= " << hnode3->Frame().GetPos().x << "    y=" << hnode3->Frame().GetPos().y << "    z=" << hnode3->Frame().GetPos().z << "\n\n";
////
////
////
////    application.DoStep();
////
////
////
////    GetLog() << "Press SPACE bar to start/stop dynamic simulation \n\n";
////    GetLog() << "Press F10 for nonlinear static solution \n\n";
////    GetLog() << "Press F11 for linear static solution \n\n";
////
////    while (application.GetDevice()->run())
////    {
////        application.BeginScene();
////
////        application.DrawAll();
////
////        application.DoStep();
////
////        application.EndScene();
////    }
////
////
////    return 0;
////}
////
////

int model_reduction_ChCSR3Matrix_OnlyRef()
{
    const int n = 10;
    chrono::ChCSR3Matrix matA(n, n);

    matA.SetElement(0, 0, 1.11);
    matA.SetElement(1, 1, 2);
    matA.SetElement(2, 2, 3);
    matA.SetElement(3, 3, 4);
    matA.SetElement(4, 4, 5);
    matA.SetElement(5, 5, 6);
    matA.SetElement(6, 6, 7);
    matA.SetElement(7, 7, 8);
    matA.SetElement(8, 8, 9);
    matA.SetElement(9, 9, 10);
    matA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChCSR3Matrix matB(n, n);

    matB.SetElement(0, 0, 1.11);
    matB.SetElement(1, 1, 2);
    matB.SetElement(2, 2, 3);
    matB.SetElement(3, 3, 4);
    matB.SetElement(4, 4, 5);
    matB.SetElement(5, 5, 6);
    matB.SetElement(6, 6, 7);
    matB.SetElement(7, 7, 8);
    matB.SetElement(8, 8, 9);
    matB.SetElement(9, 9, 10);
    matB.Compress();

    Eigen::Ref<Eigen::SparseMatrix<double>> matA_Ref = Eigen::Map<Eigen::SparseMatrix<double>>(matA.GetNumRows(), 
                                                                                                matA.GetNumColumns(),
                                                                                                matA.GetNNZ(),
                                                                                                matA.GetCSR_LeadingIndexArray(),
                                                                                                matA.GetCSR_TrailingIndexArray(),
                                                                                                matA.GetCSR_ValueArray());

    Eigen::Ref<Eigen::SparseMatrix<double>> matB_Ref = Eigen::Map<Eigen::SparseMatrix<double>>(matB.GetNumRows(),
                                                                                                matB.GetNumColumns(),
                                                                                                matB.GetNNZ(),
                                                                                                matB.GetCSR_LeadingIndexArray(),
                                                                                                matB.GetCSR_TrailingIndexArray(),
                                                                                                matB.GetCSR_ValueArray());


    //for (auto row_sel = 0; row_sel<matA_test.rows(); ++row_sel)
    //{
    //    for (auto col_sel = 0; col_sel<matA_test.cols(); ++col_sel)
    //    {
    //        std::cout << matA_test.coeff(row_sel, col_sel) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    std::cout << "************* OnlyRef test *************" << std::endl;
    std::cout << "original address: " << matA.GetCSR_ValueArray() << std::endl;
    std::cout << "OnlyRef address: " << matA_Ref.valuePtr() << std::endl;

    // Construct matrix operation object using the wrapper classes
    SparseSymMatProd<double> op(matA_Ref);
    SparseCholesky<double>  Bop(matB_Ref);
    std::cout << "coeff of OnlyRef: " << matA_Ref.coeff(0, 0) << std::endl;

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    ////// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}

int model_reduction_ChCSR3Matrix_CppRefRef()
{
    const int n = 10;
    chrono::ChCSR3Matrix matA(n, n);

    matA.SetElement(0, 0, 1.11);
    matA.SetElement(1, 1, 2);
    matA.SetElement(2, 2, 3);
    matA.SetElement(3, 3, 4);
    matA.SetElement(4, 4, 5);
    matA.SetElement(5, 5, 6);
    matA.SetElement(6, 6, 7);
    matA.SetElement(7, 7, 8);
    matA.SetElement(8, 8, 9);
    matA.SetElement(9, 9, 10);
    matA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChCSR3Matrix matB(n, n);

    matB.SetElement(0, 0, 1.11);
    matB.SetElement(1, 1, 2);
    matB.SetElement(2, 2, 3);
    matB.SetElement(3, 3, 4);
    matB.SetElement(4, 4, 5);
    matB.SetElement(5, 5, 6);
    matB.SetElement(6, 6, 7);
    matB.SetElement(7, 7, 8);
    matB.SetElement(8, 8, 9);
    matB.SetElement(9, 9, 10);
    matB.Compress();

    Eigen::Ref<Eigen::SparseMatrix<double>> matA_Ref = Eigen::Map<Eigen::SparseMatrix<double>>(matA.GetNumRows(),
                                                                                               matA.GetNumColumns(),
                                                                                               matA.GetNNZ(),
                                                                                               matA.GetCSR_LeadingIndexArray(),
                                                                                               matA.GetCSR_TrailingIndexArray(),
                                                                                               matA.GetCSR_ValueArray());

    Eigen::Ref<Eigen::SparseMatrix<double>> matB_Ref = Eigen::Map<Eigen::SparseMatrix<double>>(matB.GetNumRows(),
                                                                                               matB.GetNumColumns(),
                                                                                               matB.GetNNZ(),
                                                                                               matB.GetCSR_LeadingIndexArray(),
                                                                                               matB.GetCSR_TrailingIndexArray(),
                                                                                               matB.GetCSR_ValueArray());


    const Eigen::SparseMatrix<double>&matA_CppRef_Ref = matA_Ref;
    const Eigen::SparseMatrix<double>&matB_CppRef_Ref = matB_Ref;

    //for (auto row_sel = 0; row_sel<matA_test.rows(); ++row_sel)
    //{
    //    for (auto col_sel = 0; col_sel<matA_test.cols(); ++col_sel)
    //    {
    //        std::cout << matA_test.coeff(row_sel, col_sel) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    std::cout << "************* CppRef of Ref test *************" << std::endl;
    std::cout << "original address: " << matA.GetCSR_ValueArray() << std::endl;
    std::cout << "CppRef of Ref address: " << matA_CppRef_Ref.valuePtr() << std::endl;

    // Construct matrix operation object using the wrapper classes
    SparseSymMatProd<double> op(matA_CppRef_Ref);
    SparseCholesky<double>  Bop(matB_CppRef_Ref);
    std::cout << "coeff of CppRef of Ref: " << matA_CppRef_Ref.coeff(0, 0) << std::endl;

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    ////// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}

int model_reduction_ChCSR3Matrix_CppRefMap()
{
    const int n = 10;
    chrono::ChCSR3Matrix matA(n, n);

    matA.SetElement(0, 0, 1.11);
    matA.SetElement(1, 1, 2);
    matA.SetElement(2, 2, 3);
    matA.SetElement(3, 3, 4);
    matA.SetElement(4, 4, 5);
    matA.SetElement(5, 5, 6);
    matA.SetElement(6, 6, 7);
    matA.SetElement(7, 7, 8);
    matA.SetElement(8, 8, 9);
    matA.SetElement(9, 9, 10);
    matA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChCSR3Matrix matB(n, n);

    matB.SetElement(0, 0, 1.11);
    matB.SetElement(1, 1, 2);
    matB.SetElement(2, 2, 3);
    matB.SetElement(3, 3, 4);
    matB.SetElement(4, 4, 5);
    matB.SetElement(5, 5, 6);
    matB.SetElement(6, 6, 7);
    matB.SetElement(7, 7, 8);
    matB.SetElement(8, 8, 9);
    matB.SetElement(9, 9, 10);
    matB.Compress();

    const Eigen::SparseMatrix<double>& matA_CppRef = Eigen::Map<Eigen::SparseMatrix<double>>(matA.GetNumRows(),
                                                                                               matA.GetNumColumns(),
                                                                                               matA.GetNNZ(),
                                                                                               matA.GetCSR_LeadingIndexArray(),
                                                                                               matA.GetCSR_TrailingIndexArray(),
                                                                                               matA.GetCSR_ValueArray());

    const Eigen::SparseMatrix<double>& matB_CppRef = Eigen::Map<Eigen::SparseMatrix<double>>(matB.GetNumRows(),
                                                                                               matB.GetNumColumns(),
                                                                                               matB.GetNNZ(),
                                                                                               matB.GetCSR_LeadingIndexArray(),
                                                                                               matB.GetCSR_TrailingIndexArray(),
                                                                                               matB.GetCSR_ValueArray());


    //for (auto row_sel = 0; row_sel<matA_test.rows(); ++row_sel)
    //{
    //    for (auto col_sel = 0; col_sel<matA_test.cols(); ++col_sel)
    //    {
    //        std::cout << matA_test.coeff(row_sel, col_sel) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    std::cout << "************* CppRef test *************" << std::endl;
    std::cout << "original address: " << matA.GetCSR_ValueArray() << std::endl;
    std::cout << "CppRef of Map address: " << matA_CppRef.valuePtr() << std::endl;

    // Construct matrix operation object using the wrapper classes
    SparseSymMatProd<double> op(matA_CppRef);
    SparseCholesky<double>  Bop(matB_CppRef);
    std::cout << "coeff of CppRef: " << matA_CppRef.coeff(0, 0) << std::endl;


    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    ////// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}

int model_reduction_ChCSR3Matrix_OnlyMap()
{
    const int n = 10;
    chrono::ChCSR3Matrix matA(n, n);

    matA.SetElement(0, 0, 1.11);
    matA.SetElement(1, 1, 2);
    matA.SetElement(2, 2, 3);
    matA.SetElement(3, 3, 4);
    matA.SetElement(4, 4, 5);
    matA.SetElement(5, 5, 6);
    matA.SetElement(6, 6, 7);
    matA.SetElement(7, 7, 8);
    matA.SetElement(8, 8, 9);
    matA.SetElement(9, 9, 10);
    matA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChCSR3Matrix matB(n, n);

    matB.SetElement(0, 0, 1.11);
    matB.SetElement(1, 1, 2);
    matB.SetElement(2, 2, 3);
    matB.SetElement(3, 3, 4);
    matB.SetElement(4, 4, 5);
    matB.SetElement(5, 5, 6);
    matB.SetElement(6, 6, 7);
    matB.SetElement(7, 7, 8);
    matB.SetElement(8, 8, 9);
    matB.SetElement(9, 9, 10);
    matB.Compress();

    Eigen::Map<Eigen::SparseMatrix<double>>matA_map(matA.GetNumRows(),
                                                    matA.GetNumColumns(),
                                                    matA.GetNNZ(),
                                                    matA.GetCSR_LeadingIndexArray(),
                                                    matA.GetCSR_TrailingIndexArray(),
                                                    matA.GetCSR_ValueArray());

    Eigen::Map<Eigen::SparseMatrix<double>>matB_map(matB.GetNumRows(),
                                                    matB.GetNumColumns(),
                                                    matB.GetNNZ(),
                                                    matB.GetCSR_LeadingIndexArray(),
                                                    matB.GetCSR_TrailingIndexArray(),
                                                    matB.GetCSR_ValueArray());


    //for (auto row_sel = 0; row_sel<matA_test.rows(); ++row_sel)
    //{
    //    for (auto col_sel = 0; col_sel<matA_test.cols(); ++col_sel)
    //    {
    //        std::cout << matA_test.coeff(row_sel, col_sel) << " ";
    //    }
    //    std::cout << std::endl;
    //}
    std::cout << "************* OnlyMap test *************" << std::endl;
    std::cout << "original address: " << matA.GetCSR_ValueArray() << std::endl;
    std::cout << "OnlyMap address: " << matA_map.valuePtr() << std::endl;

    // Construct matrix operation object using the wrapper classes
    SparseSymMatProd<double> op(matA_map);
    SparseCholesky<double>  Bop(matB_map);
    std::cout << "coeff of OnlyMap: " << matA_map.coeff(0, 0) << std::endl;

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    ////// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}

int test_ChEigenSparseMatrix()
{
    chrono::ChEigenSparseMatrix mat(3, 3);

    mat.SetElement(0, 0, 1.1);
    mat.SetElement(1, 1, 2.2, false);
    mat.coeffRef(2, 2) = 3.3;

    assert(mat.rows() == mat.GetNumRows());
    assert(mat.cols() == mat.GetNumColumns());

    //for (auto row_sel = 0; row_sel<mat.rows(); ++row_sel)
    //{
    //    for (auto col_sel = 0; col_sel<mat.cols(); ++col_sel)
    //    {
    //        std::cout << mat.GetElement(2, 2, 3.3);
    //    }
    //}

    std::cout << mat.GetElement(0,0) << std::endl;
    std::cout << mat.coeff(0,0) << std::endl;



    return false;
}

int model_reduction_ChEigenSparseMatrixWrapper()
{

    const int n = 10;
    chrono::ChEigenSparseMatrixWrapper mymatA(n, n);

    mymatA.SetElement(0, 0, 1);
    mymatA.SetElement(1, 1, 2);
    mymatA.SetElement(2, 2, 3);
    mymatA.SetElement(3, 3, 4);
    mymatA.SetElement(4, 4, 5);
    mymatA.SetElement(5, 5, 6);
    mymatA.SetElement(6, 6, 7);
    mymatA.SetElement(7, 7, 8);
    mymatA.SetElement(8, 8, 9);
    mymatA.SetElement(9, 9, 10);
    mymatA.Compress();


    chrono::ChEigenSparseMatrixWrapper mymatB(n, n);

    mymatB.SetElement(0, 0, 1);
    mymatB.SetElement(1, 1, 2);
    mymatB.SetElement(2, 2, 3);
    mymatB.SetElement(3, 3, 4);
    mymatB.SetElement(4, 4, 5);
    mymatB.SetElement(5, 5, 6);
    mymatB.SetElement(6, 6, 7);
    mymatB.SetElement(7, 7, 8);
    mymatB.SetElement(8, 8, 9);
    mymatB.SetElement(9, 9, 10);
    mymatB.Compress();

    // Construct matrix operation object using the wrapper classes
    auto& mymatA_wrapper = mymatA.GetInternalEigenMatrix();
    auto& mymatB_wrapper = mymatB.GetInternalEigenMatrix();
    SparseSymMatProd<double> op(mymatA_wrapper);
    SparseCholesky<double>  Bop(mymatB_wrapper);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY>
        geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();
    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;
    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }
    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    // Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA_wrapper, mymatB_wrapper);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}

int model_reduction_ChEigenSparseMatrixWrapper_ChMatrix(chrono::ChMatrixDynamic<double>& eig_val)
{

    const int n = 10;
    chrono::ChEigenSparseMatrixWrapper mymatA(n, n);

    mymatA.SetElement(0, 0, 1);
    mymatA.SetElement(1, 1, 2);
    mymatA.SetElement(2, 2, 3);
    mymatA.SetElement(3, 3, 4);
    mymatA.SetElement(4, 4, 5);
    mymatA.SetElement(5, 5, 6);
    mymatA.SetElement(6, 6, 7);
    mymatA.SetElement(7, 7, 8);
    mymatA.SetElement(8, 8, 9);
    mymatA.SetElement(9, 9, 10);
    mymatA.Compress();


    chrono::ChEigenSparseMatrixWrapper mymatB(n, n);

    mymatB.SetElement(0, 0, 1);
    mymatB.SetElement(1, 1, 2);
    mymatB.SetElement(2, 2, 3);
    mymatB.SetElement(3, 3, 4);
    mymatB.SetElement(4, 4, 5);
    mymatB.SetElement(5, 5, 6);
    mymatB.SetElement(6, 6, 7);
    mymatB.SetElement(7, 7, 8);
    mymatB.SetElement(8, 8, 9);
    mymatB.SetElement(9, 9, 10);
    mymatB.Compress();

    // Construct matrix operation object using the wrapper classes
    auto& mymatA_wrapper = mymatA.GetInternalEigenMatrix();
    auto& mymatB_wrapper = mymatB.GetInternalEigenMatrix();
    SparseSymMatProd<double> op(mymatA_wrapper);
    SparseCholesky<double>  Bop(mymatB_wrapper);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    int num1 = 3;
    int num2 = 6;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY>
        geigs(&op, &Bop, num1, num2);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    //Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    eig_val.Resize(num1, 1);
    Eigen::Map<Eigen::VectorXd> eig_val_map(eig_val.GetAddress(), num1, 1);

    if (geigs.info() == SUCCESSFUL)
    {
        eig_val_map = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << eig_val_map << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    // Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA_wrapper, mymatB_wrapper);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}


int test_EigenMap()
{
    chrono::ChCSR3Matrix matCSR(3, 3);
    matCSR.SetElement(0, 0, 1);
    matCSR.SetElement(1, 1, 2);
    matCSR.SetElement(2, 2, 3);
    matCSR.Compress();
    std::cout << "original address: " << matCSR.GetCSR_ValueArray() << std::endl;

    Eigen::Map <Eigen::SparseMatrix<double>> eigen_matrix(matCSR.GetNumRows(), matCSR.GetNumColumns(), matCSR.GetNNZ(), matCSR.GetCSR_LeadingIndexArray(), matCSR.GetCSR_TrailingIndexArray(), matCSR.GetCSR_ValueArray());
    std::cout << "create Map: " << eigen_matrix.coeff(1, 1) << std::endl;

    Eigen::Ref<Eigen::SparseMatrix<double>> mat_Ref = Eigen::Map<Eigen::SparseMatrix<double>>(matCSR.GetNumRows(), matCSR.GetNumColumns(), matCSR.GetNNZ(), matCSR.GetCSR_LeadingIndexArray(), matCSR.GetCSR_TrailingIndexArray(), matCSR.GetCSR_ValueArray());
    std::cout << "Ref of Map: " << mat_Ref.valuePtr() << std::endl;

    const Eigen::SparseMatrix<double>& mat_pure = Eigen::Map<Eigen::SparseMatrix<double>>(matCSR.GetNumRows(), matCSR.GetNumColumns(), matCSR.GetNNZ(), matCSR.GetCSR_LeadingIndexArray(), matCSR.GetCSR_TrailingIndexArray(), matCSR.GetCSR_ValueArray());
    std::cout << "const SparseMatrix& of Map: " << mat_pure.valuePtr() << std::endl;

    const Eigen::SparseMatrix<double>& mat_Ref_ref(mat_Ref);
    std::cout << "const SparseMatrix& of Ref: " << mat_Ref_ref.valuePtr() << std::endl;

    matCSR.SetElement(1, 1, 9);
    std::cout << "coeff of const SparseMatrix& of Ref: " << mat_Ref_ref.coeff(1, 1) << std::endl;
    std::cout << "coeff of Ref of Map: " << mat_Ref.coeff(1, 1) << std::endl;



    return 0;
}


int model_reduction_ChEigenSparseMatrix()
{

    const int n = 10;
    chrono::ChEigenSparseMatrix mymatA(n, n);

    mymatA.SetElement(0, 0, 1);
    mymatA.SetElement(1, 1, 2);
    mymatA.SetElement(2, 2, 3);
    mymatA.SetElement(3, 3, 4);
    mymatA.SetElement(4, 4, 5);
    mymatA.SetElement(5, 5, 6);
    mymatA.SetElement(6, 6, 7);
    mymatA.SetElement(7, 7, 8);
    mymatA.SetElement(8, 8, 9);
    mymatA.SetElement(9, 9, 10);
    mymatA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChEigenSparseMatrix mymatB(n, n);

    mymatB.SetElement(0, 0, 1);
    mymatB.SetElement(1, 1, 2);
    mymatB.SetElement(2, 2, 3);
    mymatB.SetElement(3, 3, 4);
    mymatB.SetElement(4, 4, 5);
    mymatB.SetElement(5, 5, 6);
    mymatB.SetElement(6, 6, 7);
    mymatB.SetElement(7, 7, 8);
    mymatB.SetElement(8, 8, 9);
    mymatB.SetElement(9, 9, 10);
    mymatB.Compress();

    // Construct matrix operation object using the wrapper classes
    SparseSymMatProd<double> op(mymatA);
    SparseCholesky<double>  Bop(mymatB);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
    SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, 3, 6);

    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();

    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;

    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }

    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    ////// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = mymatB;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> es(mymatA, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}


int test_ChModelReduction()
{
    const int n = 10;
    chrono::ChCSR3Matrix matA(n, n);

    matA.SetElement(0, 0, 1.11);
    matA.SetElement(1, 1, 2);
    matA.SetElement(2, 2, 3);
    matA.SetElement(3, 3, 4);
    matA.SetElement(4, 4, 5);
    matA.SetElement(5, 5, 6);
    matA.SetElement(6, 6, 7);
    matA.SetElement(7, 7, 8);
    matA.SetElement(8, 8, 9);
    matA.SetElement(9, 9, 10);
    matA.Compress();


    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    chrono::ChCSR3Matrix matB(n, n);

    matB.SetElement(0, 0, 1.11);
    matB.SetElement(1, 1, 2);
    matB.SetElement(2, 2, 3);
    matB.SetElement(3, 3, 4);
    matB.SetElement(4, 4, 5);
    matB.SetElement(5, 5, 6);
    matB.SetElement(6, 6, 7);
    matB.SetElement(7, 7, 8);
    matB.SetElement(8, 8, 9);
    matB.SetElement(9, 9, 10);
    matB.Compress();

    chrono::ChMatrixDynamic<double> eig_val;
    chrono::ChMatrixDynamic<double> eig_vect;

    chrono::ChSymGEigsSolver eig_solver(matA, matB, eig_val, eig_vect);
    eig_solver.compute(3);

    chrono::GetLog() << eig_val << "\n";
    chrono::GetLog() << eig_vect << "\n";

    return 0;

}

int main ()
{
    //test_ChEigenSparseMatrix();
    //model_reduction_ChEigenSparseMatrixWrapper();
    //model_reduction_ChEigenSparseMatrix();
    //chrono::ChMatrixDynamic<double> eig_val;
    //model_reduction_ChEigenSparseMatrixWrapper_ChMatrix(eig_val);

    //test_EigenMap();
    int error = 0;

    try {
        model_reduction_ChCSR3Matrix_CppRefRef();
    }
    catch (...)
    {
        std::cout << "model_reduction_ChCSR3Matrix_CppRefRef broken" << std::endl;
        error += 1;
    }

    try{
        model_reduction_ChCSR3Matrix_CppRefMap();
    }
    catch(...)
    {
        std::cout << "model_reduction_ChCSR3Matrix_CppRefMap broken" << std::endl;
        error += 1;
    }

    try {
        model_reduction_ChCSR3Matrix_OnlyRef();
    }
    catch (...)
    {
        std::cout << "model_reduction_ChCSR3Matrix_OnlyRef broken" << std::endl;
        error += 1;
    }

    try {
        model_reduction_ChCSR3Matrix_OnlyMap();
    }
    catch (...)
    {
        std::cout << "model_reduction_ChCSR3Matrix_OnlyMap broken" << std::endl;
        error += 1;
    }
    


    test_ChModelReduction();



    return error;
}