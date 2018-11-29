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
// Authors: Dario Mangoni
// =============================================================================

#ifndef CHSOLVERCVXOPT_H
#define CHSOLVERCVXOPT_H

#include "chrono/core/ChTimer.h"
#include "chrono/solver/ChSolver.h"
#include "chrono/solver/ChSystemDescriptor.h"

#include "chrono_cvxopt/ChApiCvxopt.h"
#include "chrono_cvxopt/ChCvxoptConeQpEngine.h"

namespace chrono {

/// \class ChSolverCvxoptConeQp
/// Class that leverages the CVXOPT library in order to solve Chrono problems.
/// It can solve linear systems. It cannot solve VI and complementarity problems.
/// According to http://cvxopt.org/userguide/coneprog.html#quadratic-cone-programs,
/// CVXOPT solves the following problem:
/// min 0.5*x'*P*x+q'*x
/// s.t. Gx<=h
///      Ax=b

/// <pre>
///  | P   A'  G' |*|x |   |0 |   |-q|
///  | A   0   0  | |y | + |0 | = |b |
///  | G   0   0  | |z |   |s |   |h |
/// </pre>
/// with s>=0
/// This means, in particular, that unilateral constraints - removing the slack variable - are of the type Gx<=h.
/// However, Chrono uses the notation with the >= operator. We have to consider this during the conversion (swap sign of
/// h and G).
///
/// Chrono has the following structure and symbols
/// <pre>
///  | H  -C' -D' |*|q |   | 0 |   | f |
///  | C   0   0  | |le| + |-ce| = |-be|
///  | D   0  -E  | |li|   |-ci|   |-bi|
/// </pre>
///   where 'ce' is actually zero;
/// thus having this symbols equivalence
/// CVXOPT | P | A |  G | x | y | z | s | q | b | h |
/// Chrono | H | C | -D | q | le| li| ci|-f | -be| bi|

class ChApiCvxopt ChSolverCvxoptConeQp : public ChSolver {
  public:
    ChSolverCvxoptConeQp();
    virtual ~ChSolverCvxoptConeQp() {}

    /// Reset timers for internal phases in Solve and Setup.
    void ResetTimers();

    /// Enable/disable locking the sparsity pattern (default: false).\n
    /// If \a val is set to true, then the sparsity pattern of the problem matrix is assumed
    /// to be unchanged from call to call.
    void SetSparsityPatternLock(bool val);

    /// Call an update of the sparsity pattern on the underlying matrix.\n
    /// It is used to inform the solver (and the underlying matrices) that the sparsity pattern is changed.\n
    /// It is suggested to call this function just after the construction of the solver.
    /// \remark Turn on the sparsity pattern lock feature #SetSparsityPatternLock(); otherwise performance can be
    /// compromised.
    void ForceSparsityPatternUpdate(bool val = true) { m_force_sparsity_pattern_update = val; }

    /// Get cumulative time for assembly operations in Solve phase.
    double GetTimeSolve_Assembly() const { return m_timer_solve_assembly(); }
    /// Get cumulative time for Pardiso calls in Solve phase.
    double GetTimeSolve_SolverCall() const { return m_timer_solve_solvercall(); }
    /// Get cumulative time for assembly operations in Setup phase.
    double GetTimeSetup_Assembly() const { return m_timer_setup_assembly(); }
    /// Get cumulative time for Pardiso calls in Setup phase.
    double GetTimeSetup_SolverCall() const { return m_timer_setup_solvercall(); }
    /// Return the number of calls to the solver's Setup function.
    int GetNumSetupCalls() const { return m_setup_call; }
    /// Return the number of calls to the solver's Setup function.
    int GetNumSolveCalls() const { return m_solve_call; }

    bool Setup(ChSystemDescriptor& sysd) override;

    /// Solve using CVXOPT ConeQp solver.
    double Solve(ChSystemDescriptor& sysd) override;

    bool SolveRequiresMatrix() const override { return true; }

    ChCvxoptConeQpEngine& GetEngine() { return m_engine; }

  private:
    ChCvxoptConeQpEngine m_engine;  ///< interface to Cvxopt solver
    ChCSMatrixIndexType<int_t> P = {0, 0, false, 0}, G = {0, 0, false, 0}, A = {0, 0, false, 0};
    ChMatrixDynamic<double> q, h, b;     ///< rhs vector
    ChMatrixDynamic<double> x, y, z;     ///< unknowns vector
    ChMatrixDynamic<double> fric_coeff;  ///< friction coefficients

    int m_solution_dim = 0;  ///< size of the solution vector (also P rows and cols, G cols)
    unsigned int m_constr_bilat_dim = 0;
    unsigned int m_constr_unilin_dim = 0;
    unsigned int m_constr_unicone_dim = 0;
    int m_nnz_P = 0;       ///< user-supplied estimate of NNZ for matrix P
    int m_nnz_A = 0;       ///< user-supplied estimate of NNZ for matrix P
    int m_nnz_G = 0;       ///< user-supplied estimate of NNZ for matrix P
    int m_solve_call = 0;  ///< counter for calls to Solve
    int m_setup_call = 0;  ///< counter for calls to Setup

    bool m_lock = false;                           ///< is the matrix sparsity pattern locked?
    bool m_force_sparsity_pattern_update = false;  ///< is the sparsity pattern changed compared to last call?

    static const int cone_dimension = 3;

    ChTimer<> m_timer_setup_assembly;    ///< timer for matrix assembly
    ChTimer<> m_timer_setup_solvercall;  ///< timer for factorization
    ChTimer<> m_timer_solve_assembly;    ///< timer for RHS assembly
    ChTimer<> m_timer_solve_solvercall;  ///< timer for solution
};

}  // namespace chrono

#endif
