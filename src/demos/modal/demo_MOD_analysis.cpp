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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================
//
// Show how to use the ChModalAssembly to do a basic modal analysis (eigenvalues
// and eigenvector of the ChModalAssembly, which can also contain constraints.
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/core/ChRandom.h"

#include "chrono/fea/ChMesh.h"

#include "chrono/core/ChTimer.h"

#include "chrono_modal/ChModalAssembly.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "chrono/solver/ChDirectSolverLScomplex.h"
#ifdef CHRONO_PARDISO_MKL
    #include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#endif

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::modal;
using namespace chrono::fea;
using namespace chrono::irrlicht;

// Output directory
const std::string out_dir = GetChronoOutputPath() + "MODAL_ANALYSIS";

int ID_current_example = 1;

double beam_Young = 100.e6;
double beam_density = 1000;
double beam_wz = 0.3;
double beam_wy = 0.05;
double beam_L = 6;

unsigned int num_modes = 14;
bool USE_GRAVITY = true;

void MakeAndRunDemoCantilever(ChSystem& sys, ChVisualSystemIrrlicht& vis, bool base_fixed) {
    // Clear previous demo, if any:
    sys.Clear();
    sys.SetChTime(0);

    // CREATE THE ASSEMBLY.
    //
    // The ChModalAssembly is the most important item when doing modal analysis.
    // You must add finite elements, bodies and constraints into this assembly in order
    // to compute the modal frequencies etc.; objects not added into this won't be counted.
    auto assembly = chrono_types::make_shared<ChModalAssembly>();
    assembly->SetModalAutomaticGravity(USE_GRAVITY);
    sys.Add(assembly);

    // Now populate the assembly to analyze.
    // In this demo, make a cantilever constrained to a body (a large box
    // acting as a base):

    // BODY: the base:

    auto my_body_A = chrono_types::make_shared<ChBodyEasyBox>(1, 2, 2, 200);
    my_body_A->SetFixed(base_fixed);
    my_body_A->SetPos(ChVector3d(-0.5, 0, 0));
    assembly->Add(my_body_A);

    // Create a FEM mesh, that is a container for groups of elements and their referenced nodes.
    auto mesh = chrono_types::make_shared<ChMesh>();
    assembly->Add(mesh);

    mesh->SetAutomaticGravity(USE_GRAVITY);

    // BEAMS:

    // Create a simplified section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.
    auto section = chrono_types::make_shared<ChBeamSectionEulerAdvanced>();

    section->SetDensity(beam_density);
    section->SetYoungModulus(beam_Young);
    section->SetShearModulusFromPoisson(0.31);
    section->SetRayleighDampingBeta(0.00001);
    section->SetRayleighDampingAlpha(0.001);
    section->SetAsRectangularSection(beam_wy, beam_wz);

    // This helps creating sequences of nodes and ChElementBeamEuler elements:
    ChBuilderBeamEuler builder;

    builder.BuildBeam(mesh,                      // the mesh where to put the created nodes and elements
                      section,                   // the ChBeamSectionEuler to use for the ChElementBeamEuler elements
                      6,                         // the number of ChElementBeamEuler to create
                      ChVector3d(0, 0, 0),       // the 'A' point in space (beginning of beam)
                      ChVector3d(beam_L, 0, 0),  // the 'B' point in space (end of beam)
                      ChVector3d(0, 1, 0)        // the 'Y' up direction of the section for the beam
    );

    // CONSTRAINT: connect root of blade to the base.

    auto my_root = chrono_types::make_shared<ChLinkMateGeneric>();
    my_root->Initialize(builder.GetLastBeamNodes().front(), my_body_A, ChFrame<>(ChVector3d(0, 0, 1), QUNIT));
    assembly->Add(my_root);

    // VISUALIZATION ASSETS:

    auto visualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_TX);
    visualizebeamA->SetColorscaleMinMax(-0.001, 1200);
    visualizebeamA->SetSmoothFaces(true);
    visualizebeamA->SetWireframe(false);
    mesh->AddVisualShapeFEA(visualizebeamA);

    auto visualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    visualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visualizebeamC->SetSymbolsThickness(0.2);
    visualizebeamC->SetSymbolsScale(0.1);
    visualizebeamC->SetZbufferHide(false);
    mesh->AddVisualShapeFEA(visualizebeamC);

    // Just for later reference, dump M,R,K,Cq matrices. Ex. for comparison with Matlab eigs()
    // assembly->WriteSubassemblyMatrices(true, true, true, true, out_dir + "/dump");

    // Here we perform the modal analysis on the ChModalAssembly.
    // - We compute only the first n modes. This helps dealing with very large
    //   systems with many DOFs.
    // - If the assembly is free floating (ex like an helicopter or an airplane,
    //   i.e. there is no part that is fixed to ground) it will give six modes with 0 frequency,
    //   the so called rigid body modes.
    // - After computing the modes, you can access the eigenmodes, eigenvalues (also scaled as frequencies)
    //   from the ChModalAssembly member data, ex. via assembly->GetUndampedFrequencies
    // - For an interactive display of the modes, in Irrlicht view, use application.SetModalShow(true);
    //   this will pause the dynamic simulation and plot the modes of any ChModalAssembly present in the system
    //   as an oscillating animation. Use the GUI of Irrlicht 3D view to change the ID and amplitude of the plotted
    //   mode.
    ChTimer timer1;
    timer1.reset();

    timer1.start();
    ChGeneralizedEigenvalueSolverKrylovSchur eig_solver;
    ChModalSolverUndamped mod_solver(num_modes, 1e-5,
                                    500,    // max iterations per each {modes,freq} pair
                                    1e-10,  // tolerance
                                    false,  // verbose
                                    eig_solver);
    assembly->ComputeModes(mod_solver);
    timer1.stop();
    std::cout << "Time for ComputeModes: " << timer1() << std::endl;
    std::cout << "         MatrixAssembly: " << eig_solver.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "         EigenSetup: " << eig_solver.GetTimeEigenSetup() << "s" << std::endl;
    std::cout << "         GetTimeEigenSolver: " << eig_solver.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "         SolutionPostProcessing: " << eig_solver.GetTimeSolutionPostProcessing() << "s" << std::endl;

    timer1.reset();

    timer1.start();
    ChGeneralizedEigenvalueSolverKrylovSchur eig_solver2;
    ChModalSolverUndamped mod_solver2(num_modes, 1e-5,
                                     500,    // max iterations per each {modes,freq} pair
                                     1e-10,  // tolerance
                                     false,  // verbose
                                     eig_solver2);
    assembly->ComputeModesFaster(mod_solver2);
    timer1.stop();
    std::cout << "Time for ComputeModesFaster: " << timer1() << std::endl;
    std::cout << "         MatrixAssembly: " << eig_solver2.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "         EigenSetup: " << eig_solver2.GetTimeEigenSetup() << "s" << std::endl;
    std::cout << "         GetTimeEigenSolver: " << eig_solver2.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "         SolutionPostProcessing: " << eig_solver2.GetTimeSolutionPostProcessing() << "s" << std::endl;

    // Just for logging the frequencies:
    // std::cout << " The undamped modal frequencies of the modal assembly (not the whole system) in full state are: "
    //          << std::endl;
    // for (int i = 0; i < assembly->GetUndampedFrequencies().rows(); ++i) {
    //    std::cout << "Mode n." << i << "  frequency [Hz]: " << assembly->GetUndampedFrequencies()(i)
    //              << "  damping ratio:" << assembly->GetDampingRatios()(i)
    //              << "    Re=" << assembly->GetEigenValues()(i).real()
    //              << "  Im=" << assembly->GetEigenValues()(i).imag() << std::endl;
    //}

    // Here we perform the complex-modal analysis (damped modes) on the ChModalAssembly.
    // Short way:
    ////assembly->ComputeModesDamped(num_modes);
    // Or, if you need more control on the eigenvalue solver, do this:
    ////assembly->ComputeModesDamped(ChModalSolverDamped( ...parameters...));

//    // Set linear solver
//#ifdef CHRONO_PARDISO_MKL
//    ChSolverComplexPardisoMKL factorization;
//    factorization.GetMklEngine().pardisoParameterArray()[12] = 1;  // custom setting for Pardiso
//#else
//    ChSolverSparseComplexQR factorization;
//#endif
//
//    assembly->ComputeModesDamped(ChModalSolverDamped(
//        num_modes,                                              // n. of requested eigenmodes
//        1e-5,                                                   // base frequency, or vector of frequency spans
//        500,                                                    // the max. number of iterations
//        1e-10,                                                  // the tolerance
//        true,                                                   // verbose
//        ChUnsymGenEigenvalueSolverKrylovSchur(&factorization)  // the eigenpair solver algorithm
//        ));
//
//    // Just for logging the frequencies:
//    std::cout << " The damped modal results of the modal assembly (not the whole system) in full state are: "
//              << std::endl;
//    for (int i = 0; i < assembly->GetUndampedFrequencies().rows(); ++i) {
//        std::cout << "Damped mode n." << i << "  frequency [Hz]: " << assembly->GetUndampedFrequencies()(i)
//                  << "  damping ratio:" << assembly->GetDampingRatios()(i)
//                  << "    Re=" << assembly->GetEigenValues()(i).real()
//                  << "  Im=" << assembly->GetEigenValues()(i).imag() << std::endl;
//    }
//
//    // This is needed if you want to see things in Irrlicht
//    vis.BindAll();
//
//    int current_example = ID_current_example;
//    while ((ID_current_example == current_example) && vis.Run()) {
//        vis.BeginScene();
//        vis.Render();
//        tools::drawGrid(&vis, 1, 1, 12, 12, ChCoordsys<>(ChVector3d(0, 0, 0), CH_PI_2, VECT_Z),
//                        ChColor(0.5f, 0.5f, 0.5f), true);
//        vis.EndScene();
//    }
}

void MakeAndRunDemoMesh(ChSystem& sys, ChVisualSystemIrrlicht& vis) {
    // Clear previous demo, if any:
    sys.Clear();
    sys.SetChTime(0);

    // CREATE THE ASSEMBLY.
    //
    // The ChModalAssembly is the most important item when doing modal analysis.
    // You must add finite elements, bodies and constraints into this assembly in order
    // to compute the modal frequencies etc.; objects not added into this won't be counted.
    auto assembly = chrono_types::make_shared<ChModalAssembly>();
    sys.Add(assembly);

    // FINITE ELEMENT MESH
    // Create a mesh, that is a container for groups
    // of FEA elements and their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();

    // Create a material, that must be assigned to each solid element in the mesh,
    // and set its parameters
    auto mmaterial = chrono_types::make_shared<ChContinuumElastic>();
    mmaterial->SetYoungModulus(0.016e9);  // rubber 0.01e9, steel 200e9
    mmaterial->SetPoissonRatio(0.4);
    mmaterial->SetRayleighDampingBeta(0.004);
    mmaterial->SetDensity(1000);

    // Global parameter for tire:
    double tire_rad = 0.8;
    double tire_vel_z0 = -3;
    ChVector3d tire_center(0, 0.02 + tire_rad, 0.5);
    ChMatrix33<> tire_alignment(QuatFromAngleY(CH_PI));  // create rotated 180 deg on y

    // Create the surface material
    auto mysurfmaterial = chrono_types::make_shared<ChContactMaterialSMC>();
    mysurfmaterial->SetYoungModulus(10e4);
    mysurfmaterial->SetFriction(0.3f);
    mysurfmaterial->SetRestitution(0.2f);
    mysurfmaterial->SetAdhesion(0);

    // Load an ABAQUS .INP tetrahedron mesh file from disk, defining a tetrahedron mesh.
    // Note that not all features of INP files are supported. Also, quadratic tetrahedrons are promoted to linear.
    // This is much easier than creating all nodes and elements via C++ programming.
    // Ex. you can generate these .INP files using Abaqus or exporting from the SolidWorks simulation tool.

    std::map<std::string, std::vector<std::shared_ptr<ChNodeFEAbase>>> node_sets;

    // try {
    //     ChMeshFileLoader::FromAbaqusFile(my_mesh,
    //                                      GetChronoDataFile("models/tractor_wheel/tractor_wheel_coarse.INP").c_str(),
    //                                      mmaterial, node_sets, tire_center, tire_alignment);

    //} catch (std::exception myerr) {
    //    std::cerr << myerr.what() << std::endl;
    //    return;
    //}

    for (int i = 0; i < 3; ++i) {
        ChCoordsys<> cdown1(VNULL, QuatFromAngleX(CH_PI_2));
        ChCoordsys<> cdown2(ChVector3d(0, -0.4, 0));
        ChCoordsys<> crot(VNULL, QuatFromAngleY(5 * CH_2PI * ChRandom::Get()));
        for (int j = 0; j < 3; ++j) {
            try {
                ChCoordsys<> cydisp(ChVector3d(0.3 * j, 0.1 + i * 0.1, 0));
                ChCoordsys<> ctot = cdown2 >> cdown1 >> cydisp >> crot;
                ChMatrix33<> mrot(ctot.rot);
                ChMeshFileLoader::FromTetGenFile(my_mesh, GetChronoDataFile("fea/beam.node").c_str(),
                                                 GetChronoDataFile("fea/beam.ele").c_str(), mmaterial, ctot.pos, mrot);
            } catch (std::exception myerr) {
                std::cerr << myerr.what();
                return;
            }
        }
    }

    // try {
    //     ChMeshFileLoader::FromTetGenFile(my_mesh, GetChronoDataFile("fea/beam.node").c_str(),
    //                                      GetChronoDataFile("fea/beam.ele").c_str(), mmaterial);
    //     std::cout << "Mesh loaded succesfully" << std::endl;
    // } catch (std::exception myerr) {
    //     std::cerr << myerr.what() << std::endl;
    //     return;
    // }

    // Remember to add the mesh to the system!
    assembly->Add(my_mesh);

    // VISUALIZATION ASSETS:

    // auto visualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    // visualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_TX);
    // visualizebeamA->SetColorscaleMinMax(-0.001, 1200);
    // visualizebeamA->SetSmoothFaces(true);
    // visualizebeamA->SetWireframe(false);
    // mesh->AddVisualShapeFEA(visualizebeamA);

    auto visualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    visualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    visualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visualizebeamC->SetSymbolsThickness(0.2);
    visualizebeamC->SetSymbolsScale(0.1);
    visualizebeamC->SetZbufferHide(false);
    my_mesh->AddVisualShapeFEA(visualizebeamC);

    //// Just for later reference, dump M,R,K,Cq matrices. Ex. for comparison with Matlab eigs()
    // assembly->WriteSubassemblyMatrices(true, true, true, true, out_dir + "/dump");

    // Here we perform the modal analysis on the ChModalAssembly.
    // - We compute only the first n modes. This helps dealing with very large
    //   systems with many DOFs.
    // - If the assembly is free floating (ex like an helicopter or an airplane,
    //   i.e. there is no part that is fixed to ground) it will give six modes with 0 frequency,
    //   the so called rigid body modes.
    // - After computing the modes, you can access the eigenmodes, eigenvalues (also scaled as frequencies)
    //   from the ChModalAssembly member data, ex. via assembly->GetUndampedFrequencies
    // - For an interactive display of the modes, in Irrlicht view, use application.SetModalShow(true);
    //   this will pause the dynamic simulation and plot the modes of any ChModalAssembly present in the system
    //   as an oscillating animation. Use the GUI of Irrlicht 3D view to change the ID and amplitude of the plotted
    //   mode.
    ChTimer timer1;
    timer1.reset();

    timer1.start();
    ChGeneralizedEigenvalueSolverKrylovSchur eig_solver;
    ChModalSolverUndamped mod_solver(num_modes, 1e-5,
                                    500,    // max iterations per each {modes,freq} pair
                                    1e-10,  // tolerance
                                    false,  // verbose
                                    eig_solver);
    assembly->ComputeModes(mod_solver);
    timer1.stop();
    std::cout << "Time for ComputeModes: " << timer1() << std::endl;
    std::cout << "         ASS:Setup: " << assembly->GetTimeSetup() << "s" << std::endl;
    std::cout << "         ASS:MatrixAssembly: " << assembly->GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "         ASS:ModalSolver: " << assembly->GetTimeModalSolver() << "s" << std::endl;
    std::cout << "             MOD:MatrixAssembly: " << mod_solver.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "             MOD:GetTimeEigenSolver: " << mod_solver.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "                 EIG:MatrixAssembly: " << eig_solver.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "                 EIG:EigenSetup: " << eig_solver.GetTimeEigenSetup() << "s" << std::endl;
    std::cout << "                 EIG:GetTimeEigenSolver: " << eig_solver.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "                 EIG:SolutionPostProcessing: " << eig_solver.GetTimeSolutionPostProcessing() << "s"
              << std::endl;
    std::cout << "             MOD:SolutionPostProcessing: " << mod_solver.GetTimeSolutionPostProcessing() << "s"
              << std::endl;

    timer1.reset();
    assembly->ResetTimers();

    ChVectorDynamic<std::complex<double>> eigenvalues_splitmatrix(assembly->GetEigenValues());
    ChMatrixDynamic<std::complex<double>> eigenvectors_splitmatrix(assembly->GetEigenVectors());

    timer1.start();
    ChGeneralizedEigenvalueSolverKrylovSchur eig_solver2;
    ChModalSolverUndamped mod_solver2(num_modes, 1e-5,
                                     500,    // max iterations per each {modes,freq} pair
                                     1e-10,  // tolerance
                                     false,  // verbose
                                     eig_solver2);
    assembly->ComputeModesFaster(mod_solver2);
    timer1.stop();
    std::cout << "Time for ComputeModesFaster: " << timer1() << std::endl;
    std::cout << "         ASS:Setup: " << assembly->GetTimeSetup() << "s" << std::endl;
    std::cout << "         ASS:MatrixAssembly: " << assembly->GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "         ASS:ModalSolver: " << assembly->GetTimeModalSolver() << "s" << std::endl;
    std::cout << "             MOD:MatrixAssembly: " << mod_solver2.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "             MOD:GetTimeEigenSolver: " << mod_solver2.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "                 EIG:MatrixAssembly: " << eig_solver2.GetTimeMatrixAssembly() << "s" << std::endl;
    std::cout << "                 EIG:EigenSetup: " << eig_solver2.GetTimeEigenSetup() << "s" << std::endl;
    std::cout << "                 EIG:GetTimeEigenSolver: " << eig_solver2.GetTimeEigenSolver() << "s" << std::endl;
    std::cout << "                 EIG:SolutionPostProcessing: " << eig_solver2.GetTimeSolutionPostProcessing() << "s"
              << std::endl;
    std::cout << "             MOD:SolutionPostProcessing: " << mod_solver2.GetTimeSolutionPostProcessing() << "s"
              << std::endl;

    ChVectorDynamic<std::complex<double>> eigenvalues_onematrix = assembly->GetEigenValues();
    ChMatrixDynamic<std::complex<double>> eigenvectors_onematrix = assembly->GetEigenVectors();

    std::cout << "Difference eigvals: " << (eigenvalues_splitmatrix - eigenvalues_onematrix).lpNorm<Eigen::Infinity>()
              << std::endl;

    std::cout << "Difference eigvect: " << (eigenvectors_splitmatrix - eigenvectors_onematrix).lpNorm<Eigen::Infinity>()
              << std::endl;

    //    // Set linear solver
    // #ifdef CHRONO_PARDISO_MKL
    //    ChSolverComplexPardisoMKL factorization;
    //    factorization.GetMklEngine().pardisoParameterArray()[12] = 1;  // custom setting for Pardiso
    // #else
    //    ChSolverSparseComplexQR factorization;
    // #endif
    //
    //    assembly->ComputeModesDamped(ChModalSolverDamped(
    //        num_modes,                                              // n. of requested eigenmodes
    //        1e-5,                                                   // base frequency, or vector of frequency spans
    //        500,                                                    // the max. number of iterations
    //        1e-10,                                                  // the tolerance
    //        true,                                                   // verbose
    //        ChUnsymGenEigenvalueSolverKrylovSchur(&factorization)  // the eigenpair solver algorithm
    //        ));
}

void MakeAndRunDemoLbeam(ChSystem& sys, ChVisualSystemIrrlicht& vis, bool body1fixed, bool body2fixed) {
    // Clear previous demo, if any:
    sys.Clear();
    sys.SetChTime(0);

    // CREATE THE ASSEMBLY.
    // The ChModalAssembly is the most important item when doing modal analysis.
    // You must add finite elements, bodies and constraints into this assembly in order
    // to compute the modal frequencies etc.; objects not added into this won't be counted.
    auto assembly = chrono_types::make_shared<ChModalAssembly>();
    assembly->SetModalAutomaticGravity(USE_GRAVITY);
    sys.Add(assembly);

    // Now populate the assembly to analyze.
    // In this demo, make a L-shaped beam constrained to two bodies at the end

    // Create a FEM mesh, that is a container for groups of elements and their referenced nodes.
    auto mesh = chrono_types::make_shared<ChMesh>();
    assembly->Add(mesh);

    mesh->SetAutomaticGravity(USE_GRAVITY);

    // BEAMS:

    // Create a simplified section, i.e. thickness and material properties
    // for beams. This will be shared among some beams.
    auto section = chrono_types::make_shared<ChBeamSectionEulerAdvanced>();

    section->SetDensity(beam_density);
    section->SetYoungModulus(beam_Young);
    section->SetShearModulusFromPoisson(0.31);
    section->SetRayleighDampingBeta(0.00001);
    section->SetRayleighDampingAlpha(0.001);
    section->SetAsRectangularSection(beam_wy, beam_wz);

    // This helps creating sequences of nodes and ChElementBeamEuler elements:
    ChBuilderBeamEuler builder;

    builder.BuildBeam(mesh,                      // the mesh where to put the created nodes and elements
                      section,                   // the ChBeamSectionEuler to use for the ChElementBeamEuler elements
                      6,                         // the number of ChElementBeamEuler to create
                      ChVector3d(0, 0, 0),       // the 'A' point in space (beginning of beam)
                      ChVector3d(beam_L, 0, 0),  // the 'B' point in space (end of beam)
                      ChVector3d(0, 1, 0)        // the 'Y' up direction of the section for the beam
    );
    auto start_node = builder.GetLastBeamNodes().front();
    builder.BuildBeam(mesh,     // the mesh where to put the created nodes and elements
                      section,  // the ChBeamSectionEuler to use for the ChElementBeamEuler elements
                      6,        // the number of ChElementBeamEuler to create
                      builder.GetLastBeamNodes().back(),    // the 'A' point in space (beginning of beam)
                      ChVector3d(beam_L, beam_L * 0.5, 0),  // the 'B' point in space (end of beam)
                      ChVector3d(1, 0, 0)                   // the 'Y' up direction of the section for the beam
    );
    auto end_node = builder.GetLastBeamNodes().back();

    // BODY: 1st end

    auto my_body_A = chrono_types::make_shared<ChBodyEasyBox>(0.5, 0.5, 0.5, 200);
    my_body_A->SetFixed(body1fixed);
    my_body_A->SetPos(ChVector3d(-0.25, 0, 0));
    assembly->Add(my_body_A);

    // BODY: 2nd end

    auto my_body_B = chrono_types::make_shared<ChBodyEasyBox>(0.5, 0.5, 0.5, 200);
    my_body_B->SetFixed(body2fixed);
    my_body_B->SetPos(ChVector3d(beam_L, beam_L * 0.5 + 0.25, 0));
    assembly->Add(my_body_B);

    // CONSTRAINT: connect beam end to body
    auto my_root1 = chrono_types::make_shared<ChLinkMateGeneric>();
    my_root1->Initialize(start_node, my_body_A, ChFrame<>(ChVector3d(0, 0, 0), QUNIT));
    assembly->Add(my_root1);

    // CONSTRAINT: connect beam end to body
    auto my_root2 = chrono_types::make_shared<ChLinkMateGeneric>();
    my_root2->Initialize(end_node, my_body_B, ChFrame<>(ChVector3d(beam_L, beam_L * 0.5, 0), QUNIT));
    assembly->Add(my_root2);

    // VISUALIZATION ASSETS:

    auto visualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_TX);
    visualizebeamA->SetColorscaleMinMax(-0.001, 1200);
    visualizebeamA->SetSmoothFaces(true);
    visualizebeamA->SetWireframe(false);
    mesh->AddVisualShapeFEA(visualizebeamA);

    auto visualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    visualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visualizebeamC->SetSymbolsThickness(0.2);
    visualizebeamC->SetSymbolsScale(0.1);
    visualizebeamC->SetZbufferHide(false);
    mesh->AddVisualShapeFEA(visualizebeamC);

    // Just for later reference, dump M,R,K,Cq matrices. Ex. for comparison with Matlab eigs()
    assembly->WriteSubassemblyMatrices(true, true, true, true, out_dir + "/dump");

    // Here we perform the modal analysis on the ChModalAssembly.
    // - We compute only the first n modes. This helps dealing with very large
    //   systems with many DOFs.
    // - If the assembly is free floating (ex like an helicopter or an airplane,
    //   i.e. there is no part that is fixed to ground) it will give six modes with 0 frequency,
    //   the so called rigid body modes.
    // - After computing the modes, you can access the eigenmodes, eigenvalues (also scaled as frequencies)
    //   from the ChModalAssembly member data, ex. via assembly->GetUndampedFrequencies
    // - For an interactive display of the modes, in Irrlicht view, use application.SetModalShow(true);
    //   this will pause the dynamic simulation and plot the modes of any ChModalAssembly present in the system
    //   as an oscillating animation. Use the GUI of Irrlicht 3D view to change the ID and amplitude of the plotted
    //   mode.

    ChTimer timer1;
    timer1.start();
    assembly->ComputeModes(num_modes);
    timer1.stop();
    std::cout << "Time for ComputeModes: " << timer1() << std::endl;

    timer1.reset();

    timer1.start();
    assembly->ComputeModesFaster(num_modes);
    timer1.stop();
    std::cout << "Time for ComputeModesFaster: " << timer1() << std::endl;

    // If you need to enter more detailed settings for the eigenvalue solver, do this :
    /*
    assembly->ComputeModes(ChModalSolverUndamped(
        num_modes,      // n. lowest nodes to search, or modes clusters {{freq1,nnodes2},{freq2,nnodes2},{...,...}}
        1e-5,           // base freq.
        500,            // max iterations
        1e-10,          // tolerance
        false,          // verbose
        ChGeneralizedEigenvalueSolverKrylovSchur()) // solver type
    );
    */

    // Just for logging the frequencies:
    std::cout << " The modal frequencies of the modal assembly (not the whole system) in full state are: " << std::endl;
    for (int i = 0; i < assembly->GetUndampedFrequencies().rows(); ++i)
        std::cout << "Mode n." << i << "  frequency [Hz]: " << assembly->GetUndampedFrequencies()(i) << std::endl;

    // This is needed if you want to see things in Irrlicht 3D view.
    vis.BindAll();

    // Set the limitation of visualization on the mode orders.
    vis.SetModalModesMax(num_modes - 1);

    int current_example = ID_current_example;
    while ((ID_current_example == current_example) && vis.Run()) {
        vis.BeginScene();
        vis.Render();
        tools::drawGrid(&vis, 1, 1, 12, 12, ChCoordsys<>(ChVector3d(0, 0, 0), CH_PI_2, VECT_Z),
                        ChColor(0.5f, 0.5f, 0.5f), true);
        vis.EndScene();
    }
}

// Custom event manager
class MyEventReceiver : public irr::IEventReceiver {
  public:
    MyEventReceiver() {}

    bool OnEvent(const irr::SEvent& event) {
        // check if user presses keys
        if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown) {
            switch (event.KeyInput.Key) {
                case irr::KEY_KEY_1:
                    ID_current_example = 1;
                    return true;
                case irr::KEY_KEY_2:
                    ID_current_example = 2;
                    return true;
                case irr::KEY_KEY_3:
                    ID_current_example = 3;
                    return true;
                case irr::KEY_KEY_4:
                    ID_current_example = 4;
                    return true;
                case irr::KEY_KEY_5:
                    ID_current_example = 5;
                    return true;
                case irr::KEY_KEY_6:
                    ID_current_example = 6;
                    return true;
                default:
                    break;
            }
        }
        return false;
    }
};

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2021 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    // Directory for output data
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    // CREATE THE MODEL

    // Create a Chrono::Engine physical system
    ChSystemNSC sys;

    // no gravity used here
    sys.SetGravitationalAcceleration(VNULL);

    // VISUALIZATION

    // Create the Irrlicht visualization system
    ChVisualSystemIrrlicht vis;
    vis.AttachSystem(&sys);
    vis.SetWindowSize(1024, 768);
    vis.SetWindowTitle("Modal analysis");
    vis.Initialize();
    vis.AddLogo();
    vis.AddSkyBox();
    vis.AddCamera(ChVector3d(1, 1.3, 6), ChVector3d(3, 0, 0));
    vis.AddLightWithShadow(ChVector3d(20, 20, 20), ChVector3d(0, 0, 0), 50, 5, 50, 55);
    vis.AddLight(ChVector3d(-20, -20, 0), 6, ChColor(0.6f, 1.0f, 1.0f));
    vis.AddLight(ChVector3d(0, -20, -20), 6, ChColor(0.6f, 1.0f, 1.0f));

    // This is for GUI tweaking of system parameters..
    MyEventReceiver receiver;
    // note how to add a custom event receiver to the default interface:
    vis.AddUserEventReceiver(&receiver);

    // Some help on the screen
    vis.GetGUIEnvironment()->addStaticText(
        L"Press 1: fixed cantilever\n"
        L"Press 2: free-free cantilever\n"
        L"Press 3: L-beam, root fixed\n"
        L"Press 4: L-beam, free-free\n"
        L"Press 5: L-beam, fixed-fixed",
        irr::core::rect<irr::s32>(400, 80, 850, 200), false, true, 0);

    // Note, in order to have this modal visualization  working, a ChModalAssembly must have been added to the ChSystem,
    // where some modes must have been already computed.
    vis.EnableModalAnalysis(true);
    vis.SetModalSpeed(15);
    vis.SetModalAmplitude(0.8);
    vis.SetModalModeNumber(0);

    // Optional: open the GUI for changing the N of the mode shape via a slider in the Irrlicht view
    vis.ShowInfoPanel(true);
    vis.SetInfoTab(1);

    // Run the sub-demos
    // while (true) {
    // switch (ID_current_example) {
    //    case 1:
    //        MakeAndRunDemoCantilever(sys, vis,
    //                                 true);  // root fixed
    //        break;
    //    case 2:
    //        MakeAndRunDemoCantilever(sys, vis,
    //                                 false);  // root free, for free-free mode
    //        break;
    //    case 3:
    //        MakeAndRunDemoLbeam(sys, vis,
    //                            true,    // end 1 fixed
    //                            false);  // end 2 free
    //        break;
    //    case 4:
    //        MakeAndRunDemoLbeam(sys, vis,
    //                            false,   // end 1 free
    //                            false);  // end 2 free
    //        break;
    //    case 5:
    //        MakeAndRunDemoLbeam(sys, vis,
    //                            true,   // end 1 fixed
    //                            true);  // end 2 fixed
    //        break;
    //    case 6:
    MakeAndRunDemoMesh(sys, vis);  // end 2 fixed
    //        break;
    //    default:
    //        break;
    //}

    //    if (!vis.Run())
    //        break;
    //}

    return 0;
}
