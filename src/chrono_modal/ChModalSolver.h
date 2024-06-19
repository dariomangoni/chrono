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

#ifndef CHMODALSOLVER_H
#define CHMODALSOLVER_H

#include "chrono_modal/ChApiModal.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTimer.h"
#include "chrono/physics/ChAssembly.h"

#include <complex>

namespace chrono {

namespace modal {

/// Class for computing eigenvalues/eigenvectors for the DAMPED constrained system.
/// It dispatches the settings to some solver of ChUnsymGenEigenvalueSolver class.
/// It handles multiple runs of the solver if one wants to find specific ranges of frequencies.
/// Finally it guarantees that eigenvalues are sorted in the appropriate order of increasing eigenvalue modulus.
class ChApiModal ChModalSolver {
  public:
    struct ChFreqSpan {
        int nmodes;
        double freq;
    };

    /// Constructor for the case of N lower modes.
    /// Ex. ChModalSolver(7);
    /// finds first 7 lowest damped modes, using default settings (i.e. the ChUnsymGenEigenvalueSolverNullspaceDirect
    /// solver). Ex.
    ///  ChModalSolver(5, 1e-5, 500, 1e-10, false, ChUnsymGenEigenvalueSolverKrylovSchur());
    /// finds first 5 lowest damped modes using the ChUnsymGenEigenvalueSolverKrylovSchur() solver.
    ChModalSolver(int n_lower_modes,  ///< n of lower modes
                  double base_freq,   ///< frequency to whom the nodes are clustered. Use 1e-5 to get n lower modes. As
                                      ///< sigma in shift&invert, as: sigma = -pow(base_freq * CH_2PI, 2). Too small
                                      ///< gives ill conditioning (no convergence). Too large misses rigid body modes.
                  bool scaleCq,
                  bool verbose)
        : m_freq_spans({{n_lower_modes, base_freq}}),
          m_clip_position_coords(true),
          m_scaleCq(scaleCq),
          m_verbose(verbose){};

    /// Constructor for the case of multiple spans of frequency analysis
    /// ex. ChModalSolver({{10,1e-5,},{5,40}} , 500);
    /// finds first 10 lower modes, then 5 modes closest to 40 Hz, etc., using
    /// multiple runs of the solver. Closest mean that some could be higher than 40Hz,
    /// others can be lower.
    /// Another example: suppose you want the 5 lowest modes, then you also are
    /// interested in 1 high frequency mode whose frequency is already know approximately,
    /// ex. 205 Hz, then you can do ChGeneralizedEigenvalueSolverGeneric({{5,1e-5,},{1,205}}, ...).
    /// Note about overlapping ranges: if n-th run finds frequencies up to X Hz, and the (n+1)-th run finds some
    /// frequency with Y Hz where Y < X, then such Y mode(s) is discarded.
    ChModalSolver(std::vector<ChFreqSpan> freq_spans, bool scaleCq, bool verbose)
        : m_freq_spans(freq_spans), m_clip_position_coords(true), m_scaleCq(scaleCq), m_verbose(verbose){};

    virtual ~ChModalSolver(){};

    int GetNumRequestedModes() const;

    void SetClipPositionCoords(bool val) { m_clip_position_coords = val; }

    const std::vector<ChFreqSpan>& GetFrequencySpans() const { return m_freq_spans; }

    /// Get cumulative time for matrix assembly.
    double GetTimeMatrixAssembly() const { return m_timer_matrix_assembly(); }

    /// Get cumulative time eigensolver solution.
    double GetTimeEigenSolver() const { return m_timer_eigen_solver(); }

    /// Get cumulative time for post-solver solution postprocessing.
    double GetTimeSolutionPostProcessing() const { return m_timer_solution_postprocessing(); }

  protected:
    mutable ChTimer m_timer_matrix_assembly;          ///< timer for matrix assembly
    mutable ChTimer m_timer_eigen_solver;             ///< timer for eigensolver solution
    mutable ChTimer m_timer_solution_postprocessing;  ///< timer for conversion of eigensolver solution

    std::vector<ChFreqSpan> m_freq_spans;

    bool m_clip_position_coords =
        true;  ///< store only the part of each eigenvector that refers to the position coordinates
    bool m_scaleCq = true;
    bool m_verbose = false;
};

}  // end namespace modal

}  // end namespace chrono

#endif
