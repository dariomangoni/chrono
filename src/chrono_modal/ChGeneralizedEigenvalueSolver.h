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

#ifndef CHGENERALIZEDEIGENVALUESOLVER_H
#define CHGENERALIZEDEIGENVALUESOLVER_H

#include "chrono_modal/ChApiModal.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTimer.h"
#include "chrono/physics/ChAssembly.h"

#include <complex>
#include <functional>
#include <numeric>

namespace chrono {

void ChApiModal
PlaceMatrix(ChSparseMatrix& HCQ, const ChSparseMatrix& H, int row_start, int col_start, double scale = 1.0);

Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor, int>> ChApiModal
getColMajorSparseMatrix(const ChSparseMatrix& mat);

void ChApiModal CountNonZerosForEachRow(const ChSparseMatrix& Q, Eigen::VectorXi& nonZerosPerRow, int offset);

void ChApiModal CountNonZerosForEachRowTransposed(const ChSparseMatrix& Q_transp,
                                                  Eigen::VectorXi& nonZerosPerRow,
                                                  int offset);

namespace modal {

template <typename EigSolverType>
int Solve(EigSolverType& eig_solver,
          ChSparseMatrix& A,
          ChSparseMatrix& B,
          ChMatrixDynamic<typename EigSolverType::ScalarType>& eigvects,
          ChVectorDynamic<typename EigSolverType::ScalarType>& eigvals,
          const std::list<std::pair<int, typename EigSolverType::ScalarType>>& eig_requests);

/// Base interface class for generalized eigenvalue solvers.
/// Currently it implies that the derived eigensolvers are iterative.
template <typename ScalarT>
class ChGeneralizedEigenvalueSolver {
  public:
    using ScalarType = typename ScalarT;

    ChGeneralizedEigenvalueSolver() {}

    virtual ~ChGeneralizedEigenvalueSolver(){};

    double tolerance = 1e-6;              ///< tolerance for the iterative solver.
    int max_iterations = 500;             ///< upper limit for the number of iterations. If too low might not converge.
    bool verbose = false;                 ///< turn to true to see some diagnostic.
    mutable bool sort_ritz_pairs = true;  ///< sort the eigenvalues based on the smallest absolute value

    /// Get cumulative time for matrix assembly.
    double GetTimeMatrixAssembly() const { return m_timer_matrix_assembly(); }

    /// Get cumulative time eigensolver setup.
    double GetTimeEigenSetup() const { return m_timer_eigen_setup(); }

    /// Get cumulative time eigensolver solution.
    double GetTimeEigenSolver() const { return m_timer_eigen_solver(); }

    /// Get cumulative time for post-solver solution postprocessing.
    double GetTimeSolutionPostProcessing() const { return m_timer_solution_postprocessing(); }

    /// Sort the eigenvalues and eigenvectors in the order specified by the ordering function in-place.
    static void SortRitzPairs(
        ChVectorDynamic<ScalarType>& eigvals,
        ChMatrixDynamic<ScalarType>& eigvects,
        std::function<bool(const ChVectorDynamic<ScalarType>&, int, int)> ordering_fun =
            [](const ChVectorDynamic<std::complex<double>>& eigv, int a, int b) {
                return std::abs(eigv(a)) < std::abs(eigv(b));
            }) {
        // TODO: do the permutation only if not already ordered
        auto perm = GetPermutationMatrix(
            eigvals.rows(), std::bind(ordering_fun, eigvals, std::placeholders::_1, std::placeholders::_2));

        eigvects *= perm;
        eigvals = perm.transpose() * eigvals;
    }

    static double GetMaxResidual(const ChSparseMatrix& A,
                                 const ChSparseMatrix& B,
                                 const ChMatrixDynamic<ScalarType>& eigvects,
                                 const ChVectorDynamic<ScalarType>& eigvals) {
        double max_residual = -1;

        ChVectorDynamic<ScalarType> residuals_col;
        for (int i = 0; i < eigvects.cols(); i++) {
            residuals_col = A * eigvects.col(i) - eigvals(i) * B * eigvects.col(i);
            double cur_residual = residuals_col.lpNorm<Eigen::Infinity>();
            max_residual = std::max(cur_residual, max_residual);
        }

        return max_residual;
    }

    static double GetMaxResidual(const ChSparseMatrix& K,
                                 const ChSparseMatrix& M,
                                 const ChSparseMatrix& Cq,
                                 const ChMatrixDynamic<ScalarType>& eigvects,
                                 const ChVectorDynamic<ScalarType>& eigvals) {
        int n = K.rows();
        int m = Cq.rows();
        double max_residual = 0;
        for (auto nv = 0; nv < eigvals.size(); nv++) {
            double cur_residual_state =
                (K * eigvects.col(nv).topRows(n) + Cq.transpose() * eigvects.col(nv).bottomRows(m) +
                 eigvals(nv) * M * eigvects.col(nv).topRows(n))
                    .lpNorm<Eigen::Infinity>();
            double cur_residual_lambda = (Cq * eigvects.col(nv).topRows(n)).lpNorm<Eigen::Infinity>();
            double cur_residual = std::max(cur_residual_state, cur_residual_lambda);
            if (cur_residual > max_residual) {
                max_residual = cur_residual;
            }
        }

        return max_residual;
    }

    static double GetMaxResidual(const ChSparseMatrix& K,
                                 const ChSparseMatrix& R,
                                 const ChSparseMatrix& M,
                                 const ChSparseMatrix& Cq,
                                 const ChMatrixDynamic<ScalarType>& eigvects,
                                 const ChVectorDynamic<ScalarType>& eigvals) {
        int n = K.rows();
        int m = Cq.rows();
        double max_residual = 0;
        for (auto nv = 0; nv < eigvals.size(); nv++) {
            double cur_residual_state_p =
                (eigvects.col(nv).segment(n, n) - eigvals(nv) * eigvects.col(nv).topRows(n)).lpNorm<Eigen::Infinity>();

            double cur_residual_state_v =
                (K * eigvects.col(nv).topRows(n) + R * eigvects.col(nv).segment(n, n) +
                 Cq.transpose() * eigvects.col(nv).bottomRows(m) + eigvals(nv) * M * eigvects.col(nv).segment(n, n))
                    .lpNorm<Eigen::Infinity>();
            double cur_residual_lambda = (Cq * eigvects.col(nv).topRows(n)).lpNorm<Eigen::Infinity>();
            double cur_residual = std::max(std::max(cur_residual_state_p, cur_residual_state_v), cur_residual_lambda);
            if (cur_residual > max_residual) {
                max_residual = cur_residual;
            }
        }

        return max_residual;
    }

    /// Build the quadratic eigenvalue problem matrix for the undamped case.
    /// The optional scaling of the Cq matrix improves the conditioning of the eigenvalue problem.
    /// If enabled, the user must take care of rescaling back the eigenvectors;
    /// e.g. eigvects.bottomRows(Cq.rows()) *= scaling
    /// The scaling value is returned by the function.
    static double BuildUndampedSystem(const ChSparseMatrix& M,
                                      const ChSparseMatrix& K,
                                      const ChSparseMatrix& Cq,
                                      ChSparseMatrix& A,
                                      ChSparseMatrix& B,
                                      bool scaleCq) {
        int n_vars = M.rows();
        int n_constr = Cq.rows();

        // Scale constraints matrix
        double scaling = 1.0;
        if (scaleCq) {
            scaling = K.diagonal().mean();
        }

        // Preallocate the A matrix
        A.resize(n_vars + n_constr, n_vars + n_constr);
        Eigen::VectorXi A_resSize;
        A_resSize.resize(A.rows());
        A_resSize.setZero();

        CountNonZerosForEachRow(K, A_resSize, 0);
        CountNonZerosForEachRowTransposed(Cq, A_resSize, 0);
        CountNonZerosForEachRow(Cq, A_resSize, n_vars);
        A.reserve(A_resSize);

        // A  =  [ -K   -Cq' ]
        //       [ -Cq    0  ]
        A.setZero();
        PlaceMatrix(A, K, 0, 0, -1.0);
        PlaceMatrix(A, Cq.transpose(), 0, n_vars, -scaling);
        PlaceMatrix(A, Cq, n_vars, 0, -scaling);
        A.makeCompressed();

        // Preallocate the B matrix
        B.resize(n_vars + n_constr, n_vars + n_constr);
        Eigen::VectorXi B_resSize;
        B_resSize.resize(B.rows());
        B_resSize.setZero();
        CountNonZerosForEachRow(M, B_resSize, 0);
        B.reserve(B_resSize);

        // B  =  [  M     0  ]
        //       [  0     0  ]

        B.setZero();
        PlaceMatrix(B, M, 0, 0);
        B.makeCompressed();

        return scaling;
    }

    /// Build the quadratic eigenvalue problem matrix for the damped case.
    /// The optional scaling of the Cq matrix improves the conditioning of the eigenvalue problem.
    /// If enabled, the user must take care of rescaling back the eigenvectors;
    /// e.g. eigvects.bottomRows(Cq.rows()) *= scaling
    /// The scaling value is returned by the function.
    static double BuildDampedSystem(const ChSparseMatrix& M,
                                    const ChSparseMatrix& R,
                                    const ChSparseMatrix& K,
                                    const ChSparseMatrix& Cq,
                                    ChSparseMatrix& A,
                                    ChSparseMatrix& B,
                                    bool scaleCq) {
        // Generate the A and B in state space
        int n_vars = M.rows();
        int n_constr = Cq.rows();

        // Scale constraints matrix
        double scaling = 1.0;
        if (scaleCq) {
            // std::cout << "Scaling Cq" << std::endl;
            scaling = K.diagonal().mean();
        }

        A.resize(2 * n_vars + n_constr, 2 * n_vars + n_constr);
        B.resize(2 * n_vars + n_constr, 2 * n_vars + n_constr);
        ChSparseMatrix identity_n_vars(n_vars, n_vars);
        identity_n_vars.setIdentity();

        // Preallocate the A matrix
        Eigen::VectorXi A_nnz_rows;
        A_nnz_rows.resize(A.rows());
        A_nnz_rows.setZero();

        A_nnz_rows.segment(0, n_vars) = ChVectorDynamic<int>::Ones(n_vars);

        CountNonZerosForEachRow(K, A_nnz_rows, n_vars);
        CountNonZerosForEachRow(R, A_nnz_rows, n_vars);
        CountNonZerosForEachRow(Cq, A_nnz_rows, 2 * n_vars);
        CountNonZerosForEachRowTransposed(Cq, A_nnz_rows, n_vars);

        A.reserve(A_nnz_rows);

        // A  =  [  0     I     0 ]
        //       [ -K    -R  -Cq' ]
        //       [ -Cq    0     0 ]

        A.setZero();
        PlaceMatrix(A, -K, n_vars, 0);
        PlaceMatrix(A, -scaling * Cq, 2 * n_vars, 0);
        PlaceMatrix(A, identity_n_vars, 0, n_vars);
        PlaceMatrix(A, -R, n_vars, n_vars);
        PlaceMatrix(A, -scaling * Cq.transpose(), n_vars, 2 * n_vars);
        A.makeCompressed();

        // Preallocate the B matrix
        Eigen::VectorXi B_nnz_rows;
        B_nnz_rows.resize(B.rows());
        B_nnz_rows.setZero();
        B_nnz_rows.segment(0, n_vars) = ChVectorDynamic<int>::Ones(n_vars);
        CountNonZerosForEachRow(M, B_nnz_rows, n_vars);

        B.reserve(B_nnz_rows);

        // B  =  [  I     0     0 ]
        //       [  0     M     0 ]
        //       [  0     0     0 ]

        B.setZero();
        PlaceMatrix(B, identity_n_vars, 0, 0);
        PlaceMatrix(B, M, n_vars, n_vars);
        B.makeCompressed();

        return scaling;
    }

    static void InsertUniqueRitzPairs(const ChVectorDynamic<ScalarType>& eigvals_source,
                                      const ChMatrixDynamic<ScalarType>& eigvects_source,
                                      ChVectorDynamic<ScalarType>& eigvals_total,
                                      ChMatrixDynamic<ScalarType>& eigvects_total,
                                      std::function<double(ScalarType)> freq_from_eigval_fun,
                                      int& found_eigs,
                                      int num_eigvals_source,
                                      double threshold = 1e-6) {
        for (auto eigv = 0; eigv < num_eigvals_source; ++eigv) {
            // look if another eigenvalue with same frequency has already been found
            bool eig_is_unique = true;
            for (auto stored_eigv_idx = 0; stored_eigv_idx < found_eigs; stored_eigv_idx++) {
                if (std::abs(freq_from_eigval_fun(eigvals_total[stored_eigv_idx]) -
                             freq_from_eigval_fun(eigvals_source[eigv])) < threshold) {
                    eig_is_unique = false;
                    break;
                }
            }
            // if the eigenvalue is unique, store it
            if (eig_is_unique) {
                eigvals_total[found_eigs] = eigvals_source[eigv];
                eigvects_total.col(found_eigs) = eigvects_source.col(eigv);
                found_eigs++;
            }
        }
    }

  protected:
    static Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> GetPermutationMatrix(
        int num_elements,
        std::function<bool(int, int)> ordering_fun) {
        std::vector<int> order(num_elements);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), ordering_fun);
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
        perm.indices() = Eigen::Map<Eigen::ArrayXi>(order.data(), order.size());
        return perm;
    }

    mutable ChTimer m_timer_matrix_assembly;          ///< timer for matrix assembly
    mutable ChTimer m_timer_eigen_setup;              ///< timer for eigensolver setup
    mutable ChTimer m_timer_eigen_solver;             ///< timer for eigensolver solution
    mutable ChTimer m_timer_solution_postprocessing;  ///< timer for conversion of eigensolver solution

    int m_min_subspace_size = 30;  // TODO: make it constant or tunable

    template <typename EigSolverType>
    friend int Solve(EigSolverType& eig_solver,
                     ChSparseMatrix& A,
                     ChSparseMatrix& B,
                     ChMatrixDynamic<typename EigSolverType::ScalarType>& eigvects,
                     ChVectorDynamic<typename EigSolverType::ScalarType>& eigvals,
                     const std::list<std::pair<int, typename EigSolverType::ScalarType>>& eig_requests);
};

template <typename EigSolverType>
int Solve(EigSolverType& eig_solver,
          ChSparseMatrix& A,
          ChSparseMatrix& B,
          ChMatrixDynamic<typename EigSolverType::ScalarType>& eigvects,
          ChVectorDynamic<typename EigSolverType::ScalarType>& eigvals,
          const std::list<std::pair<int, typename EigSolverType::ScalarType>>& eig_requests) {
    int num_modes_total = 0;
    for (const auto& eig_req : eig_requests) {
        num_modes_total += eig_req.first;
    }

    eigvals.resize(num_modes_total);
    int found_eigs = 0;

    // disable the automatic sorting of eigenvalues for any case
    bool sort_ritz_pairs_bkp = eig_solver.sort_ritz_pairs;
    eig_solver.sort_ritz_pairs = false;

    if (eig_requests.size() == 0)
        return 0;
    else if (eig_requests.size() == 1) {
        eigvects.resize(A.rows(), num_modes_total);

        eig_solver.m_timer_eigen_solver.start();

        int converged_eigs = eig_solver.Solve(A, B, eigvects, eigvals, num_modes_total, eig_requests.begin()->second);

        found_eigs = std::min(num_modes_total, converged_eigs);

        eig_solver.m_timer_eigen_solver.stop();

        eig_solver.m_timer_solution_postprocessing.start();

        if (found_eigs != num_modes_total)
            eigvects.conservativeResize(Eigen::NoChange_t::NoChange, found_eigs);

    } else {
        eigvects.resize(A.rows(), num_modes_total);

        // for each freq_spans finds the closest modes to i-th input frequency:
        for (const auto& eig_req : eig_requests) {
            eig_solver.m_timer_matrix_assembly.start();

            ChMatrixDynamic<typename EigSolverType::ScalarType> eigvects_singlespan;
            ChVectorDynamic<typename EigSolverType::ScalarType> eigvals_singlespan;

            eig_solver.m_timer_matrix_assembly.stop();

            eig_solver.m_timer_eigen_solver.start();

            int converged_eigs =
                eig_solver.Solve(A, B, eigvects_singlespan, eigvals_singlespan, eig_req.first, eig_req.second);
            eig_solver.m_timer_eigen_solver.stop();

            eig_solver.m_timer_solution_postprocessing.start();

            int min_converged_eigs = std::min(eig_req.first, converged_eigs);

            eig_solver.InsertUniqueRitzPairs(eigvals_singlespan, eigvects_singlespan, eigvals, eigvects,
                                             eig_solver.GetNaturalFrequency, found_eigs, min_converged_eigs);

            eig_solver.m_timer_solution_postprocessing.stop();
        }
        eig_solver.m_timer_solution_postprocessing.start();

        if (found_eigs != num_modes_total)
            eigvects.conservativeResize(Eigen::NoChange_t::NoChange, found_eigs);
    }

    if (found_eigs != num_modes_total)
        eigvals.conservativeResize(found_eigs);

    // re-enable original ritz pairs sorting
    eig_solver.sort_ritz_pairs = sort_ritz_pairs_bkp;

    if (eig_solver.sort_ritz_pairs)
        eig_solver.SortRitzPairs(eigvals, eigvects);

    eig_solver.m_timer_solution_postprocessing.stop();

    return found_eigs;
}

}  // end namespace modal

}  // end namespace chrono

#endif
