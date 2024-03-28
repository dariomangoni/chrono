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
// Author: Dario Mangoni
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBuilderBeam.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/timestepper/ChAssemblyAnalysis.h"

#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"

#include "chrono/fea/ChMesh.h"

#include "chrono/core/ChTimer.h"

#include "chrono_modal/ChUnsymGenEigenvalueSolver.h"
#include "chrono_modal/ChSymGenEigenvalueSolver.h"
#include "chrono_modal/ChGeneralizedEigenvalueSolver.h"


#include "chrono/solver/ChDirectSolverLScomplex.h"

#include "chrono_thirdparty/filesystem/path.h"
#include <iomanip>

// #include <unsupported/Eigen/SparseExtra>

#include <fast_matrix_market/app/Eigen.hpp>

#include "gtest/gtest.h"

using namespace chrono;
using namespace chrono::modal;
using namespace chrono::fea;

static const std::string val_dir = "../RESULTS/";
static const std::string out_dir = val_dir + "modal/";
static const std::string ref_dir = "testing/modal/analysis/";

static const double tolerance = 1e-3;

double GetEigenvaluesMaxDiff(const ChVectorDynamic<double>& eig1, const ChVectorDynamic<double>& eig2) {
    return (eig1 - eig2).lpNorm<Eigen::Infinity>();
}

double GetEigenvaluesMaxDiff(const ChVectorDynamic<std::complex<double>>& eig1,
                             const ChVectorDynamic<std::complex<double>>& eig2) {
    return std::max((eig1.real() - eig2.real()).lpNorm<Eigen::Infinity>(),
                    (eig1.imag().cwiseAbs() - eig2.imag().cwiseAbs()).lpNorm<Eigen::Infinity>());
}

void prepare_folders(std::string testname) {
    // Create output directory (if it does not already exist)
    if (!filesystem::create_directory(filesystem::path(val_dir))) {
        std::cerr << "Error creating directory " << val_dir << std::endl;
        throw std::invalid_argument("Error creating directory " + val_dir);
    }
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cerr << "Error creating directory " << out_dir << std::endl;
        throw std::invalid_argument("Error creating directory " + out_dir);
    }
    if (!filesystem::create_directory(filesystem::path(out_dir + testname))) {
        std::cerr << "Error creating directory " << out_dir + testname << std::endl;
        throw std::invalid_argument("Error creating directory " + out_dir + testname);
    }
}

std::shared_ptr<ChAssembly> BuildBeamFixBody(ChSystem& sys) {
    /*
     * Beam with end body
     *
     *   (fixed node)----()----()----()----()<--link-->[body]
     *
     */

    auto assembly = chrono_types::make_shared<ChAssembly>();

    sys.Add(assembly);

    double beam_Young = 100.e6;
    double beam_density = 1000;
    double beam_wz = 0.3;
    double beam_wy = 0.05;
    double beam_L = 6;

    double body_end_xwidth = 0.5;

    // beam
    auto mesh = chrono_types::make_shared<ChMesh>();
    assembly->Add(mesh);

    mesh->SetAutomaticGravity(false);

    auto section = chrono_types::make_shared<ChBeamSectionEulerAdvanced>();

    section->SetDensity(beam_density);
    section->SetYoungModulus(beam_Young);
    section->SetShearModulusFromPoisson(0.31);
    section->SetRayleighDampingBeta(0.00001);
    section->SetRayleighDampingAlpha(0.001);
    section->SetAsRectangularSection(beam_wy, beam_wz);

    ChBuilderBeamEuler builder;

    builder.BuildBeam(mesh,                      // the mesh where to put the created nodes and elements
                      section,                   // the ChBeamSectionEuler to use for the ChElementBeamEuler elements
                      4,                         // the number of ChElementBeamEuler to create
                      ChVector3d(0, 0, 0),       // the 'A' point in space (beginning of beam)
                      ChVector3d(beam_L, 0, 0),  // the 'B' point in space (end of beam)
                      ChVector3d(0, 1, 0)        // the 'Y' up direction of the section for the beam
    );

    builder.GetLastBeamNodes().front()->SetFixed(true);

    auto body_end = chrono_types::make_shared<ChBodyEasyBox>(body_end_xwidth, 1, 1, 200);
    body_end->SetPos(ChVector3d(beam_L + body_end_xwidth / 2.0, 0, 0));
    assembly->Add(body_end);

    auto link_beamend_body = chrono_types::make_shared<ChLinkMateFix>();
    link_beamend_body->Initialize(builder.GetLastBeamNodes().back(), body_end,
                                  ChFrame<>(ChVector3d(beam_L, 0, 0), QUNIT));
    assembly->Add(link_beamend_body);

    sys.Setup();
    sys.Update();

    return assembly;
}

void generateMRKCqfromAssembly(std::shared_ptr<ChAssembly> assembly, std::string refname) {
    ChSystemNSC sys;
    // auto assembly = BuildBeamFixBody(sys);

    ChSparseMatrix M;
    ChSparseMatrix K;
    ChSparseMatrix R;
    ChSparseMatrix Cq;

    assembly->Setup();
    assembly->Update();

    ChSystemDescriptor temp_descriptor;

    temp_descriptor.BeginInsertion();
    assembly->InjectVariables(temp_descriptor);
    assembly->InjectKRMMatrices(temp_descriptor);
    assembly->InjectConstraints(temp_descriptor);
    temp_descriptor.EndInsertion();

    temp_descriptor.UpdateCountsAndOffsets();

    // Generate the A and B in state space
    int n_vars = temp_descriptor.CountActiveVariables();
    int n_constr = temp_descriptor.CountActiveConstraints();
    M.resize(n_vars, n_vars);
    K.resize(n_vars, n_vars);
    R.resize(n_vars, n_vars);
    Cq.resize(n_constr, n_vars);

    // Stiffness matrix
    assembly->LoadKRMMatrices(-1.0, 0.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(K, 0, 0);

    // Damping matrix
    assembly->LoadKRMMatrices(0.0, 1.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(R, 0, 0);

    // Mass matrix
    assembly->LoadKRMMatrices(0.0, 0.0, 1.0);
    temp_descriptor.SetMassFactor(1.0);
    temp_descriptor.PasteMassKRMMatrixInto(M, 0, 0);

    // Constraint Jacobian
    assembly->LoadConstraintJacobians();
    temp_descriptor.PasteConstraintsJacobianMatrixInto(Cq, 0, 0);

    std::ofstream stream_M(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_M.txt"));
    std::ofstream stream_K(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_K.txt"));
    std::ofstream stream_R(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_R.txt"));
    std::ofstream stream_Cq(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Cq.txt"));

    if (stream_M.fail() || stream_K.fail() || stream_R.fail() || stream_Cq.fail()) {
        std::cerr << "Error opening file for writing in " << ref_dir + refname << " folder" << std::endl;
        return;
    }

    fast_matrix_market::write_matrix_market_eigen(stream_M, M);
    fast_matrix_market::write_matrix_market_eigen(stream_K, K);
    fast_matrix_market::write_matrix_market_eigen(stream_R, R);
    fast_matrix_market::write_matrix_market_eigen(stream_Cq, Cq);
}

// int main() {
//     std::string testname = "SymMKCqChrono";
//
//     // Create a system
//     ChSystemNSC sys;
//     auto assembly = BuildBeamFixBody(sys);
//     generateMRKCqfromAssembly(assembly, testname);
//
//     return 0;
// }

TEST(CountNonZerosForEachRow, Count) {
    std::string refname = "CountNonZeros";

    ChSparseMatrix Q;
    Eigen::VectorXi Q_nnz_rows_MATLAB;
    Eigen::VectorXi Q_nnz_cols_MATLAB;

    std::ifstream stream_Q(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Q.txt"));
    std::ifstream stream_Q_nnz_rows(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Q_nnz_rows.txt"));
    std::ifstream stream_Q_nnz_cols(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Q_nnz_cols.txt"));
    fast_matrix_market::read_matrix_market_eigen(stream_Q, Q);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_Q_nnz_rows, Q_nnz_rows_MATLAB);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_Q_nnz_cols, Q_nnz_cols_MATLAB);

    Eigen::VectorXi Q_nnz_rows_CHRONO(Q.rows());
    Q_nnz_rows_CHRONO.setZero();
    Eigen::VectorXi Q_nnz_cols_CHRONO(Q.cols());
    Q_nnz_cols_CHRONO.setZero();

    CountNonZerosForEachRow(Q, Q_nnz_rows_CHRONO, 0);
    CountNonZerosForEachRowTransposed(Q, Q_nnz_cols_CHRONO, 0);

    ASSERT_EQ(Q_nnz_rows_CHRONO, Q_nnz_rows_MATLAB);
    ASSERT_EQ(Q_nnz_cols_CHRONO, Q_nnz_cols_MATLAB);
}

template <typename EigenSolverType, typename ScalarType>
void ExecuteEigenSolverCallAB(EigenSolverType eigen_solver, std::string refname) {
    ChSparseMatrix A;
    ChSparseMatrix B;
    ChMatrixDynamic<ScalarType> sigma_mat;
    ChMatrixDynamic<int> reqeigs_mat;
    ChMatrixDynamic<ScalarType> eigvects_MATLAB;
    ChVectorDynamic<ScalarType> eigvals_MATLAB;

    ChMatrixDynamic<ScalarType> eigvects_CHRONO;
    ChVectorDynamic<ScalarType> eigvals_CHRONO;

    std::ifstream stream_A(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_A.txt"));
    std::ifstream stream_B(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_B.txt"));
    std::ifstream stream_sigma(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_sigma.txt"));
    std::ifstream stream_reqeigs(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_reqeigs.txt"));
    std::ifstream stream_eigvals_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvals_MATLAB.txt"));
    std::ifstream stream_eigvects_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvects_MATLAB.txt"));

    fast_matrix_market::read_matrix_market_eigen(stream_A, A);
    fast_matrix_market::read_matrix_market_eigen(stream_B, B);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_sigma, sigma_mat);
    ScalarType sigma = sigma_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_reqeigs, reqeigs_mat);
    int reqeigs = reqeigs_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvals_MATLAB, eigvals_MATLAB);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvects_MATLAB, eigvects_MATLAB);

    eigen_solver.Solve(A, B, eigvects_CHRONO, eigvals_CHRONO, reqeigs, sigma);

    eigen_solver.SortRitzPairs(eigvals_MATLAB, eigvects_MATLAB);

    // std::ofstream stream_eigvals_CHRONO(
    //     utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvals_CHRONO.txt"));
    // std::ofstream stream_eigvects_CHRONO(
    //     utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvects_CHRONO.txt"));
    // fast_matrix_market::write_matrix_market_eigen_dense(stream_eigvals_CHRONO, eigvals_CHRONO);
    // fast_matrix_market::write_matrix_market_eigen_dense(stream_eigvects_CHRONO, eigvects_CHRONO);

    double max_delta_eigvals = GetEigenvaluesMaxDiff(eigvals_CHRONO, eigvals_MATLAB);
    ASSERT_NEAR(max_delta_eigvals, 0, tolerance) << "Eigenvalues not matching.\n"
                                                 << "MATLAB:\n"
                                                 << eigvals_MATLAB << "\nCHRONO:\n"
                                                 << eigvals_CHRONO << std::endl;

    double max_residual_CHRONO = 0;
    int max_residual_CHRONO_idx = -1;
    for (auto nv = 0; nv < eigvals_CHRONO.size(); nv++) {
        double cur_residual =
            (A * eigvects_CHRONO.col(nv) - eigvals_CHRONO(nv) * B * eigvects_CHRONO.col(nv)).lpNorm<Eigen::Infinity>();
        if (cur_residual > max_residual_CHRONO) {
            max_residual_CHRONO = cur_residual;
            max_residual_CHRONO_idx = nv;
        }
    }

    ASSERT_NEAR(max_residual_CHRONO, 0, tolerance)
        << "Residuals exceeding threshold (index: " << max_residual_CHRONO_idx << ")" << std::endl;
}

TEST(ChSymGenEigenvalueSolverKrylovSchur, SymAB) {
    ExecuteEigenSolverCallAB<ChSymGenEigenvalueSolverKrylovSchur, double>(ChSymGenEigenvalueSolverKrylovSchur(),
                                                                          "SymAB");
}

TEST(ChSymGenEigenvalueSolverLanczos, SymAB) {
    ExecuteEigenSolverCallAB<ChSymGenEigenvalueSolverLanczos, double>(ChSymGenEigenvalueSolverLanczos(), "SymAB");
}

TEST(ChSymGenEigenvalueSolverKrylovSchur, SymMKCqChrono_AB) {
    ExecuteEigenSolverCallAB<ChSymGenEigenvalueSolverKrylovSchur, double>(ChSymGenEigenvalueSolverKrylovSchur(),
                                                                          "SymMKCqChrono");
}

TEST(ChSymGenEigenvalueSolverLanczos, SymMKCqChrono_AB) {
    ExecuteEigenSolverCallAB<ChSymGenEigenvalueSolverLanczos, double>(ChSymGenEigenvalueSolverLanczos(),
                                                                      "SymMKCqChrono");
}

void ExecuteEigenSolverCallMKCq(ChSymGenEigenvalueSolver& eigen_solver, std::string refname) {
    ChSparseMatrix M;
    ChSparseMatrix K;
    ChSparseMatrix Cq;
    ChMatrixDynamic<double> sigma_mat;
    ChMatrixDynamic<int> reqeigs_mat;
    ChMatrixDynamic<double> eigvects_MATLAB;
    ChVectorDynamic<double> eigvals_MATLAB;

    ChMatrixDynamic<double> eigvects_CHRONO;
    ChVectorDynamic<double> eigvals_CHRONO;

    std::ifstream stream_M(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_M.txt"));
    std::ifstream stream_K(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_K.txt"));
    std::ifstream stream_Cq(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Cq.txt"));
    std::ifstream stream_sigma(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_sigma.txt"));
    std::ifstream stream_reqeigs(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_reqeigs.txt"));
    std::ifstream stream_eigvals_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvals_MATLAB.txt"));
    std::ifstream stream_eigvects_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvects_MATLAB.txt"));

    fast_matrix_market::read_matrix_market_eigen(stream_M, M);
    fast_matrix_market::read_matrix_market_eigen(stream_K, K);
    fast_matrix_market::read_matrix_market_eigen(stream_Cq, Cq);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_sigma, sigma_mat);
    double sigma = sigma_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_reqeigs, reqeigs_mat);
    int reqeigs = reqeigs_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvals_MATLAB, eigvals_MATLAB);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvects_MATLAB, eigvects_MATLAB);

    eigen_solver.SortRitzPairs(eigvals_MATLAB, eigvects_MATLAB);

    // ChSymGenEigenvalueSolverKrylovSchur eigen_solver;
    const bool scaleCq = true;

    ChSparseMatrix A, B;
    ChGeneralizedEigenvalueSolver<double>::BuildUndampedSystem(M, K, Cq, A, B, scaleCq);

    eigen_solver.Solve(A, B, eigvects_CHRONO, eigvals_CHRONO, reqeigs, sigma);

    // instead of doing a simple difference, consider the imaginary part to be the same if positive or negative
    double max_delta_eigvals = GetEigenvaluesMaxDiff(eigvals_CHRONO, eigvals_MATLAB);

    ASSERT_NEAR(max_delta_eigvals, 0, tolerance) << "Eigenvalues not matching.\n"
                                                 << "MATLAB:\n"
                                                 << eigvals_MATLAB << "\nCHRONO:\n"
                                                 << eigvals_CHRONO << std::endl;
}

TEST(ChSymGenEigenvalueSolverKrylovSchur, SymMKCqChrono) {
    ChSymGenEigenvalueSolverKrylovSchur eigen_solver;
    ExecuteEigenSolverCallMKCq(eigen_solver, "SymMKCq");
}

TEST(ChSymGenEigenvalueSolverLanczos, SymMKCqChrono) {
    ChSymGenEigenvalueSolverLanczos eigen_solver;
    ExecuteEigenSolverCallMKCq(eigen_solver, "SymMKCqChrono");
}

TEST(ChSymGenEigenvalueSolverKrylovSchur, SymMKCq) {
    ChSymGenEigenvalueSolverKrylovSchur eigen_solver;
    ExecuteEigenSolverCallMKCq(eigen_solver, "SymMKCq");
}

TEST(ChSymGenEigenvalueSolverLanczos, SymMKCq) {
    ChSymGenEigenvalueSolverLanczos eigen_solver;
    ExecuteEigenSolverCallMKCq(eigen_solver, "SymMKCq");
}

TEST(ChUnsymGenEigenvalueSolverKrylovSchur, UnsymAB) {
    ExecuteEigenSolverCallAB<ChUnsymGenEigenvalueSolverKrylovSchur, std::complex<double>>(
        ChUnsymGenEigenvalueSolverKrylovSchur(chrono_types::make_shared<ChSolverSparseComplexLU>()), "UnsymAB");
}

void ExecuteEigenSolverCallMRKCq(std::string refname) {
    ChSparseMatrix M;
    ChSparseMatrix K;
    ChSparseMatrix R;
    ChSparseMatrix Cq;
    ChMatrixDynamic<std::complex<double>> sigma_mat;
    ChMatrixDynamic<int> reqeigs_mat;
    ChMatrixDynamic<std::complex<double>> eigvects_MATLAB;
    ChVectorDynamic<std::complex<double>> eigvals_MATLAB;

    ChMatrixDynamic<std::complex<double>> eigvects_CHRONO;
    ChVectorDynamic<std::complex<double>> eigvals_CHRONO;

    std::ifstream stream_M(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_M.txt"));
    std::ifstream stream_R(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_R.txt"));
    std::ifstream stream_K(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_K.txt"));
    std::ifstream stream_Cq(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Cq.txt"));
    std::ifstream stream_sigma(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_sigma.txt"));
    std::ifstream stream_reqeigs(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_reqeigs.txt"));
    std::ifstream stream_eigvals_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvals_MATLAB.txt"));
    std::ifstream stream_eigvects_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvects_MATLAB.txt"));

    fast_matrix_market::read_matrix_market_eigen(stream_M, M);
    fast_matrix_market::read_matrix_market_eigen(stream_K, K);
    fast_matrix_market::read_matrix_market_eigen(stream_R, R);
    fast_matrix_market::read_matrix_market_eigen(stream_Cq, Cq);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_sigma, sigma_mat);
    std::complex<double> sigma = sigma_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_reqeigs, reqeigs_mat);
    int reqeigs = reqeigs_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvals_MATLAB, eigvals_MATLAB);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvects_MATLAB, eigvects_MATLAB);

    ChUnsymGenEigenvalueSolverKrylovSchur eigen_solver(chrono_types::make_shared<ChSolverSparseComplexLU>());

    const bool scaleCq = true;
    ChSparseMatrix A, B;
    eigen_solver.BuildDampedSystem(M, R, K, Cq, A, B, scaleCq);

    eigen_solver.Solve(A, B, eigvects_CHRONO, eigvals_CHRONO, reqeigs, sigma);

    eigen_solver.SortRitzPairs(eigvals_MATLAB, eigvects_MATLAB);

    // instead of doing a simple difference, consider the imaginary part to be the same if positive or negative
    double max_delta_eigvals = GetEigenvaluesMaxDiff(eigvals_CHRONO, eigvals_MATLAB);

    ASSERT_NEAR(max_delta_eigvals, 0, tolerance) << "Eigenvalues not matching.\n"
                                                 << "MATLAB:\n"
                                                 << eigvals_MATLAB << "\nCHRONO:\n"
                                                 << eigvals_CHRONO << std::endl;
}

TEST(ChUnsymGenEigenvalueSolverKrylovSchur, UnsymMRKCq) {
    ExecuteEigenSolverCallMRKCq("UnsymMRKCq");
}

TEST(ChUnsymGenEigenvalueSolverKrylovSchur, UnsymMRKCq_multifreq) {
    std::string refname = "UnsymMRKCq_multifreq";

    ChSparseMatrix M;
    ChSparseMatrix K;
    ChSparseMatrix R;
    ChSparseMatrix Cq;
    ChMatrixDynamic<std::complex<double>> sigma_1_mat;
    ChMatrixDynamic<int> reqeigs_1_mat;
    ChMatrixDynamic<std::complex<double>> sigma_2_mat;
    ChMatrixDynamic<int> reqeigs_2_mat;
    ChMatrixDynamic<std::complex<double>> eigvects_MATLAB;
    ChVectorDynamic<std::complex<double>> eigvals_MATLAB;

    ChMatrixDynamic<std::complex<double>> eigvects_CHRONO;
    ChVectorDynamic<std::complex<double>> eigvals_CHRONO;

    std::ifstream stream_M(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_M.txt"));
    std::ifstream stream_R(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_R.txt"));
    std::ifstream stream_K(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_K.txt"));
    std::ifstream stream_Cq(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_Cq.txt"));
    std::ifstream stream_sigma_1(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_sigma_1.txt"));
    std::ifstream stream_sigma_2(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_sigma_2.txt"));
    std::ifstream stream_reqeigs_1(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_reqeigs_1.txt"));
    std::ifstream stream_reqeigs_2(utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_reqeigs_2.txt"));
    std::ifstream stream_eigvals_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvals_MATLAB.txt"));
    std::ifstream stream_eigvects_MATLAB(
        utils::GetValidationDataFile(ref_dir + refname + "/" + refname + "_eigvects_MATLAB.txt"));

    fast_matrix_market::read_matrix_market_eigen(stream_M, M);
    fast_matrix_market::read_matrix_market_eigen(stream_K, K);
    fast_matrix_market::read_matrix_market_eigen(stream_R, R);
    fast_matrix_market::read_matrix_market_eigen(stream_Cq, Cq);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_sigma_1, sigma_1_mat);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_sigma_2, sigma_2_mat);
    std::complex<double> sigma_1 = sigma_1_mat(0, 0);
    std::complex<double> sigma_2 = sigma_2_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_reqeigs_1, reqeigs_1_mat);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_reqeigs_2, reqeigs_2_mat);
    int reqeigs_1 = reqeigs_1_mat(0, 0);
    int reqeigs_2 = reqeigs_2_mat(0, 0);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvals_MATLAB, eigvals_MATLAB);
    fast_matrix_market::read_matrix_market_eigen_dense(stream_eigvects_MATLAB, eigvects_MATLAB);

    ChUnsymGenEigenvalueSolverKrylovSchur eigen_solver(chrono_types::make_shared<ChSolverSparseComplexLU>());

    const bool scaleCq = true;
    int eigvects_clipping_length = M.rows();

    ChSparseMatrix A, B;
    eigen_solver.BuildDampedSystem(M, R, K, Cq, A, B, scaleCq);

    std::list<std::pair<int, std::complex<double>>> eig_requests;
    eig_requests.push_back(std::make_pair(reqeigs_1, sigma_1));
    eig_requests.push_back(std::make_pair(reqeigs_2, sigma_2));
    eigen_solver.sort_ritz_pairs = true;
    modal::Solve<>(eigen_solver, A, B, eigvects_CHRONO, eigvals_CHRONO, eig_requests, eigvects_clipping_length);

    // instead of doing a simple difference, consider the imaginary part to be the same if positive or negative
    double max_delta_eigvals = GetEigenvaluesMaxDiff(eigvals_CHRONO, eigvals_MATLAB);

    ASSERT_NEAR(max_delta_eigvals, 0, tolerance) << "Eigenvalues not matching.\n"
                                                 << "MATLAB:\n"
                                                 << eigvals_MATLAB << "\nCHRONO:\n"
                                                 << eigvals_CHRONO << std::endl;
}
