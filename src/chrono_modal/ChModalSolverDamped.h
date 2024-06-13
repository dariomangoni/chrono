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

#ifndef CHMODALSOLVERDAMPED_H
#define CHMODALSOLVERDAMPED_H

#include "chrono_modal/ChApiModal.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTimer.h"
#include "chrono/physics/ChAssembly.h"
#include "chrono_modal/ChModalSolver.h"
#include "chrono_modal/ChUnsymGenEigenvalueSolver.h"

#include <complex>

namespace chrono {

// Forward references
class ChUnsymGenEigenvalueSolver;

namespace modal {

/// Class for computing eigenvalues/eigenvectors for the DAMPED constrained system.
/// It dispatches the settings to some solver of ChUnsymGenEigenvalueSolver class.
/// It handles multiple runs of the solver if one wants to find specific ranges of frequencies.
/// Finally it guarantees that eigenvalues are sorted in the appropriate order of increasing eigenvalue modulus.
/// (lambda^2*M + lambda*R + K)*x = 0 s.t. Cq*x = 0
/// also (-w^2*M + i*w*R + K)*x = 0,  with complex w (where w.length() = undamped nat.freq)
class ChApiModal ChModalSolverDamped : public ChModalSolver {
  public:
    /// Constructor for the case of N lower modes.
    /// Ex.
    ///  ChModalSolverDamped(7);
    /// finds first 7 lowest damped modes, using default settings (i.e. the ChUnsymGenEigenvalueSolverNullspaceDirect
    /// solver). Ex.
    ///  ChModalSolverDamped(5, 1e-5, 500, 1e-10, false, ChUnsymGenEigenvalueSolverKrylovSchur());
    /// finds first 5 lowest damped modes using the ChUnsymGenEigenvalueSolverKrylovSchur() solver.
    ChModalSolverDamped(
        int n_lower_modes,  ///< n of lower modes
        double base_freq,   ///< frequency to whom the nodes are clustered. Use 1e-5 to get n lower modes. As
                            ///< sigma in shift&invert, as: sigma = -pow(base_freq * CH_2PI, 2). Too small gives
                            ///< ill conditioning (no convergence). Too large misses rigid body modes.
        bool clip_position_coords = true,
        bool scaleCq = true,   ///< apply scaling to the Cq matrix to improve conditioning
        bool verbose = false,  ///< turn to true to see some diagnostic
        std::shared_ptr<ChUnsymGenEigenvalueSolver> solver = chrono_types::make_shared<
            ChUnsymGenEigenvalueSolverKrylovSchur>()  /// solver to use (default: direct, null-space based)
        )
        : ChModalSolver(n_lower_modes, base_freq, clip_position_coords, scaleCq, verbose), m_solver(solver){};

    /// Constructor for the case of multiple spans of frequency analysis
    /// ex. ChModalSolverDamped({{10,1e-5,},{5,40}} , 500);
    /// finds first 10 lower modes, then 5 modes closest to 40 Hz, etc., using
    /// multiple runs of the solver. Closest mean that some could be higher than 40Hz,
    /// others can be lower.
    /// Another example: suppose you want the 5 lowest modes, then you also are
    /// interested in 1 high frequency mode whose frequency is already know approximately,
    /// ex. 205 Hz, then you can do ChGeneralizedEigenvalueSolverGeneric({{5,1e-5,},{1,205}}, ...).
    /// Note about overlapping ranges: if n-th run finds frequencies up to X Hz, and the (n+1)-th run finds some
    /// frequency with Y Hz where Y < X, then such Y mode(s) is discarded.
    ChModalSolverDamped(
        std::vector<ChFreqSpan> freq_spans,  ///< vector of {nmodes,freq}_i , will provide first nmodes_i starting at
                                              ///< freq_i per each i vector entry
        bool clip_position_coords = true,
        bool scaleCq = true,   ///< apply scaling to the Cq matrix to improve conditioning
        bool verbose = false,  ///< turn to true to see some diagnostic
        std::shared_ptr<ChUnsymGenEigenvalueSolver> solver = chrono_types::make_shared<
            ChUnsymGenEigenvalueSolverKrylovSchur>()  /// solver to use (default: direct, null-space based)
        )
        : ChModalSolver(freq_spans, clip_position_coords, scaleCq, verbose), m_solver(solver){};

    virtual ~ChModalSolverDamped(){};

    /// Solve the constrained eigenvalue problem (-wsquare*M + K)*x = 0 s.t. Cq*x = 0
    /// Return the number of found modes, where n is not necessarily n_lower_modes (or the sum of ChFreqSpan::nmodes if
    /// multiple spans)
    virtual int Solve(
        const ChAssembly& assembly,  ///< input M matrix, n_v x n_v
        ChMatrixDynamic<std::complex<double>>&
            eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<std::complex<double>>& eigvals,  ///< output vector with n eigenvalues, will be resized.
        ChVectorDynamic<double>& freq,  ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
        ChVectorDynamic<double>& damp_ratios  ///< output vector with n damping ratios, will be resized.
    ) const;

  protected:
    std::shared_ptr<ChUnsymGenEigenvalueSolver> m_solver;
};

}  // end namespace modal

}  // end namespace chrono

#endif
