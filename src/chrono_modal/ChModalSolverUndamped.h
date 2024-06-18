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

#ifndef CHMODALSOLVERUNDAMPED_H
#define CHMODALSOLVERUNDAMPED_H

#include "chrono_modal/ChApiModal.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTimer.h"
#include "chrono/physics/ChAssembly.h"
#include "chrono_modal/ChGeneralizedEigenvalueSolver.h"
#include "chrono_modal/ChSymGenEigenvalueSolver.h"
#include "chrono_modal/ChModalSolver.h"

namespace chrono {

namespace modal {

void BuildGeneralizedEigenProblemMatrices(ChAssembly& assembly,
                                          ChSystemDescriptor& temp_descriptor,
                                          ChSparseMatrix& A,
                                          ChSparseMatrix& B,
                                          int n_vars);

/// Class for computing eigenvalues/eigenvectors for the undamped constrained system.
/// It dispatches the settings to some solver of ChGeneralizedEigenvalueSolver class.
/// It handles multiple runs of the solver if one wants to find specific ranges of frequencies.
/// Finally it guarantees that eigenvalues are sorted in the appropriate order of increasing frequency.
template <typename EigenvalueSolverType, typename ScalarType>
class ChModalSolverUndamped : public ChModalSolver {
  public:
    /// Constructor for the case of N lower modes.
    /// Ex.
    ///  ChModalSolverUndamped(7);
    /// finds first 7 lowest modes, using default settings (i.e. the ChGeneralizedEigenvalueSolverLanczos).
    /// Ex.
    ///  ChModalSolverUndamped(5, 1e-5, 500, 1e-10, false, ChSymGenEigenvalueSolverLanczos());
    /// finds first 5 lowest modes using the ChSymGenEigenvalueSolverLanczos() solver.
    ChModalSolverUndamped(
        int n_lower_modes,        ///< n of lower modes
        double base_freq = 1e-5,  ///< frequency to whom the nodes are clustered. Use 1e-5 to get n lower modes. As
                                  ///< sigma in shift&invert, as: sigma = -pow(base_freq * CH_2PI, 2). Too small gives
                                  ///< ill conditioning (no convergence). Too large misses rigid body modes
        bool clip_position_coords = true,
        bool scaleCq = true,   ///< apply scaling to the Cq matrix to improve conditioning
        bool verbose = false,  ///< turn to true to see some diagnostic
        const EigenvalueSolverType& solver = ChSymGenEigenvalueSolverLanczos())
        : ChModalSolver(n_lower_modes, base_freq, clip_position_coords, scaleCq, verbose), m_solver(solver){};

    /// Constructor for the case of multiple spans of frequency analysis
    /// ex. ChModalSolverUndamped({{10,1e-5,},{5,40}} , 500) ;
    /// finds first 10 lower modes, then 5 modes closest to 40 Hz, etc., using
    /// multiple runs of the solver. Closest mean that some could be higher than 40Hz,
    /// others can be lower.
    /// Another example: suppose you want the 5 lowest modes, then you also are
    /// interested in 1 high frequency mode whose frequency is already know approximately,
    /// ex. 205 Hz, then you can do ChModalSolverUndamped({{5,1e-5,},{1,205}}, ...).
    /// Note about overlapping ranges: if n-th run finds frequencies up to X Hz, and the (n+1)-th run finds some
    /// frequency with Y Hz where Y < X, then such Y mode(s) is discarded.
    ChModalSolverUndamped(std::vector<ChFreqSpan> freq_spans,  ///< vector of {nmodes,freq}_i , will provide first
                                                               ///< nmodes_i starting at freq_i per each i vector entry
                          bool clip_position_coords = true,
                          bool scaleCq = true,   ///< apply scaling to the Cq matrix to improve conditioning
                          bool verbose = false,  ///< turn to true to see some diagnostic
                          const EigenvalueSolverType& solver = ChSymGenEigenvalueSolverLanczos())
        : ChModalSolver(freq_spans, clip_position_coords, scaleCq, verbose), m_solver(solver){};

    virtual ~ChModalSolverUndamped(){};

    /// Solve the constrained eigenvalue problem (-wsquare*M + K)*x = 0 s.t. Cq*x = 0
    /// Return the n. of found modes, where n is not necessarily n_lower_modes (or the sum of ChFreqSpan::nmodes if
    /// multiple spans)
    int Solve(
        const ChAssembly& assembly,             ///< assembly on which to apply the eigen solver
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with n eigenvalues, will be resized.
        ChVectorDynamic<ScalarType>& freq  ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ) const;

    int Solve(
        const ChSparseMatrix& K,
        const ChSparseMatrix& M,
        const ChSparseMatrix& Cq,
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with n eigenvalues, will be resized.
        ChVectorDynamic<ScalarType>& freq  ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ) const;

  protected:
    const EigenvalueSolverType& m_solver;
};

template <typename EigenvalueSolverType, typename ScalarType>
int ChModalSolverUndamped<EigenvalueSolverType, ScalarType>::Solve(const ChAssembly& assembly,
                                                                   ChMatrixDynamic<ScalarType>& eigvects,
                                                                   ChVectorDynamic<ScalarType>& eigvals,
                                                                   ChVectorDynamic<ScalarType>& freq) const {
    ChAssembly* assembly_nonconst = const_cast<ChAssembly*>(&assembly);

    std::shared_ptr<ChSystemDescriptor> descriptor_bkp;
    if (assembly_nonconst->GetSystem()->GetSystemDescriptor())
        descriptor_bkp = assembly_nonconst->GetSystem()->GetSystemDescriptor();

    m_timer_matrix_assembly.start();

    ChSystemDescriptor sysd;
    ChSystemDescriptor temp_descriptor;

    temp_descriptor.BeginInsertion();
    assembly_nonconst->InjectVariables(temp_descriptor);
    assembly_nonconst->InjectKRMMatrices(temp_descriptor);
    assembly_nonconst->InjectConstraints(temp_descriptor);
    temp_descriptor.EndInsertion();

    temp_descriptor.UpdateCountsAndOffsets();

    // Generate the A and B in state space
    int n_vars = temp_descriptor.CountActiveVariables();
    int n_constr = temp_descriptor.CountActiveConstraints();

    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]

    // B  =  [  M     0  ]
    //       [  0     0  ]

    ChSparsityPatternLearner A_spl(n_vars + n_constr, n_vars + n_constr);
    ChSparsityPatternLearner B_spl(n_vars + n_constr, n_vars + n_constr);
    BuildGeneralizedEigenProblemMatrices(*assembly_nonconst, temp_descriptor, A_spl, B_spl, n_vars);

    ChSparseMatrix A(n_vars + n_constr, n_vars + n_constr);
    ChSparseMatrix B(n_vars + n_constr, n_vars + n_constr);
    A_spl.Apply(A);
    B_spl.Apply(B);

    A.setZeroValues();
    B.setZeroValues();
    BuildGeneralizedEigenProblemMatrices(*assembly_nonconst, temp_descriptor, A, B, n_vars);

    m_timer_matrix_assembly.stop();

    m_timer_matrix_assembly.start();

    // Scale constraints matrix
    double scaling = 1.0;
    if (m_scaleCq) {
        scaling = 0.0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (ChSparseMatrix::InnerIterator it(A, k); it; ++it) {
                if (it.row() < n_vars && it.col() == it.row()) {
                    scaling += it.valueRef();
                }
            }
        }
        scaling = scaling / n_vars;
    }

    // TODO: check scaling!
    // Cq scaling
    for (auto row_i = n_vars; row_i < n_vars + n_constr; row_i++) {
        for (auto nnz_i = A.outerIndexPtr()[row_i];
             nnz_i <
             (A.isCompressed() ? A.outerIndexPtr()[row_i + 1] : A.outerIndexPtr()[row_i] + A.innerNonZeroPtr()[row_i]);
             ++nnz_i) {
            A.valuePtr()[nnz_i] *= -scaling;
        }
    }

    // CqT scaling
    for (unsigned int k = 0; k < n_vars; ++k) {
        for (ChSparseMatrix::InnerIterator it(A, k); it; ++it) {
            if (it.col() >= n_vars) {
                it.valueRef() *= -scaling;
            }
        }
    }

    A.makeCompressed();
    B.makeCompressed();

    std::list<std::pair<int, double>> eig_requests;
    for (int i = 0; i < m_freq_spans.size(); i++) {
        eig_requests.push_back(std::make_pair(m_freq_spans[i].nmodes, m_solver.GetOptimalShift(m_freq_spans[i].freq)));
    }

    m_timer_matrix_assembly.stop();

    m_timer_eigen_solver.start();
    int found_eigs =
        modal::Solve<>(m_solver, A, B, eigvects, eigvals, eig_requests, m_clip_position_coords ? n_vars : A.rows());
    m_timer_eigen_solver.stop();

    m_timer_solution_postprocessing.start();
    m_solver.GetNaturalFrequencies(eigvals, freq);
    m_timer_solution_postprocessing.stop();

    if (descriptor_bkp) {
        assembly_nonconst->GetSystem()->SetSystemDescriptor(descriptor_bkp);
        assembly_nonconst->Setup();
    }

    return found_eigs;
}

template <typename EigenvalueSolverType, typename ScalarType>
int ChModalSolverUndamped<EigenvalueSolverType, ScalarType>::Solve(const ChSparseMatrix& K,
                                                                   const ChSparseMatrix& M,
                                                                   const ChSparseMatrix& Cq,
                                                                   ChMatrixDynamic<ScalarType>& eigvects,
                                                                   ChVectorDynamic<ScalarType>& eigvals,
                                                                   ChVectorDynamic<ScalarType>& freq) const {
    m_timer_matrix_assembly.start();

    // Generate the A and B in state space
    int n_vars = K.rows();
    int n_constr = Cq.rows();

    ChSparseMatrix A(n_vars + n_constr, n_vars + n_constr);
    ChSparseMatrix B(n_vars + n_constr, n_vars + n_constr);
    m_solver::BuildUndampedSystem(K, M, Cq, A, B, m_scaleCq);

    std::list<std::pair<int, double>> eig_requests;
    for (int i = 0; i < m_freq_spans.size(); i++) {
        eig_requests.push_back(std::make_pair(m_freq_spans[i].nmodes, m_solver.GetOptimalShift(m_freq_spans[i].freq)));
    }

    m_timer_matrix_assembly.stop();

    m_timer_eigen_solver.start();
    int found_eigs =
        modal::Solve<>(m_solver, A, B, eigvects, eigvals, eig_requests, m_clip_position_coords ? n_vars : A.rows());
    m_timer_eigen_solver.stop();

    m_timer_solution_postprocessing.start();
    m_solver.GetNaturalFrequencies(eigvals, freq);
    m_timer_solution_postprocessing.stop();

    if (descriptor_bkp) {
        assembly_nonconst->GetSystem()->SetSystemDescriptor(descriptor_bkp);
        assembly_nonconst->Setup();
    }

    return found_eigs;
}

}  // end namespace modal

}  // end namespace chrono

#endif
