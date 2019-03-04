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

#include "chrono_cvxopt/ChSolverCvxoptConeQp.h"

namespace chrono {

bool ChSolverCvxoptConeQp::Setup(ChSystemDescriptor& sysd) {
    m_setup_call++;
    return true;
}

double ChSolverCvxoptConeQp::Solve(ChSystemDescriptor& sysd) {
    // TODO: required because at EndInsertion only the standard constraints counter is updated
    sysd.UpdateCountsAndOffsetsSplit(false);

    m_timer_solve_assembly.start();
    // Calculate problem size at first call.
    m_solution_dim = sysd.CountActiveVariables();
    m_constr_bilat_dim = sysd.GetNumConstraintsBilateral();
    m_constr_unilin_dim = sysd.GetNumConstraintsUnilateralLinear();
    m_constr_unicone_dim = sysd.GetNumConstraintsUnilateralConic();
    auto m_constr_uni_dim = m_constr_unilin_dim + m_constr_unicone_dim;

    m_engine.SetNumConstrUnilateralLinear(m_constr_unilin_dim);
    assert(m_constr_unicone_dim % cone_dimension == 0 &&
           "Friction constraints dimensions are not multiple of the given 'cone_dimension'");
    m_engine.SetNumConstrUnilateralConic(
        std::vector<Py_ssize_t>(m_constr_unicone_dim / cone_dimension, cone_dimension));

    // Let the matrix acquire the information about ChSystem
    if (m_force_sparsity_pattern_update) {
        m_force_sparsity_pattern_update = false;

        ChSparsityPatternLearner sl_P(m_solution_dim, m_solution_dim, false);
        ChSparsityPatternLearner sl_G(m_constr_uni_dim, m_solution_dim, false);
        ChSparsityPatternLearner sl_A(m_constr_bilat_dim, m_solution_dim, false);
        sysd.ConvertToMatrixForm(&sl_A, &sl_G, &sl_P, nullptr, nullptr, nullptr, nullptr, nullptr, false);
        A.LoadSparsityPattern(sl_A);
        G.LoadSparsityPattern(sl_G);
        P.LoadSparsityPattern(sl_P);
    }else {
        // If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
        // a call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
        // Otherwise, do this only at the first call, using the default sparsity fill-in.
        if (m_nnz_P != 0) {
            P.Reset(m_solution_dim, m_solution_dim, m_nnz_P);
        } else if (m_setup_call == 0) {
            P.Reset(m_solution_dim, m_solution_dim,
                    static_cast<int>(m_solution_dim * (m_solution_dim * SPM_DEF_FULLNESS)));
        }

        if (m_nnz_G != 0) {
            G.Reset(m_constr_uni_dim, m_solution_dim, m_nnz_G);
        } else if (m_setup_call == 0) {
            G.Reset(m_constr_uni_dim, m_solution_dim,
                    static_cast<int>(m_constr_uni_dim * (m_solution_dim * SPM_DEF_FULLNESS)));
        }

        if (m_nnz_A != 0) {
            A.Reset(m_constr_bilat_dim, m_solution_dim, m_nnz_A);
        } else if (m_setup_call == 0) {
            A.Reset(m_constr_bilat_dim, m_solution_dim,
                    static_cast<int>(m_constr_bilat_dim * (m_solution_dim * SPM_DEF_FULLNESS)));
        }
    }

    sysd.ConvertToMatrixForm(&A, &G, &P, nullptr, &q, &b, &h, &fric_coeff, false);

    // Modify fric_coeff so that it has the friction coefficient on U and V directions instead of N
    for (auto i = m_constr_unilin_dim; i<m_constr_uni_dim; i = i+3)
    {
        auto fric_coeff_temp = fric_coeff(i, 0);
        fric_coeff(i, 0) = 1.0;
        fric_coeff(i+1, 0) = fric_coeff_temp;
        fric_coeff(i+2, 0) = fric_coeff_temp;
    }

    // Change signs
    q.MatrScale(-1);
    b.MatrScale(-1);
    G *= -1;

    // Scale G in order to project constraints on self-dual cones
    //G.ExportToDatFile("C:/orig", 6);
    std::function<void(int, int, double*)> func = [this](int r, int c, double* v) { *v = *v * this->fric_coeff(r, 0); };
    G.ForEachExistentValueInRange(func, 0, G.GetNumRows() - 1, 0, G.GetNumColumns() - 1);
    //G.ExportToDatFile("C:/mod", 6);

    h.MatrScale(fric_coeff);

    // Allow the matrix to be compressed.
    m_timer_solve_assembly.stop();
    P.Compress();
    G.Compress();
    A.Compress();
    m_timer_solve_assembly.start();

    // Set current matrix in the MKL engine.
    m_engine.SetPMatrix(P);
    m_engine.SetGMatrix(G);
    m_engine.SetAMatrix(A);

    m_engine.SetqMatrix(q);
    m_engine.SethMatrix(h);
    m_engine.SetbMatrix(b);

    m_timer_solve_assembly.stop();

    // Solve
    m_timer_solve_solvercall.start();
    m_engine.Run();
    m_timer_solve_solvercall.stop();

    m_solve_call++;

    m_timer_solve_assembly.start();
    m_engine.GetSolution(x);
    m_engine.GetLagrangianMultBilateral(y);
    m_engine.GetLagrangianMultUnilateral(z);
    m_timer_solve_assembly.stop();

    if (!m_engine.IsSolutionOptimal())
        GetLog() << "\nCVXOPT: solution at solver call: " << m_solve_call << " is not optimal.\n";

    if (verbose) {
        GetLog() << " CVXOPT Solve call: " << m_solve_call << "\n";
        GetLog() << "  Assembly: " << m_timer_solve_assembly.GetTimeSecondsIntermediate() << "s"
                 << "  Solve: " << m_timer_solve_solvercall.GetTimeSecondsIntermediate() << "s\n";
    }

    // Scale lagrangian multipliers
    z.MatrScale(fric_coeff);
    sysd.SetDynamicVariables(x);
    sysd.SetConstrReactionVariables(y, z);

    //ChMatrixDynamic<double> s;
    //m_engine.GetSlackVariable(s);
    //s.MatrDivScale(fric_coeff);

    // GetLog() << "G" << "\n";
    // PrintMatrix(G);
    // GetLog() << "\n";
    // GetLog() << fric_coeff << "\n";
    // GetLog() << "z" << z << "\n";
    // GetLog() << "s" << s << "\n";

    // system("cls");

    return 0.0;
}

void ChSolverCvxoptConeQp::SetSparsityPatternLock(bool val) {
    m_lock = val;
    P.SetSparsityPatternLock(m_lock);
    G.SetSparsityPatternLock(m_lock);
    A.SetSparsityPatternLock(m_lock);
}

ChSolverCvxoptConeQp::ChSolverCvxoptConeQp() {
    m_engine.SetShowProgress(false);
}

void ChSolverCvxoptConeQp::ResetTimers() {
    m_timer_setup_assembly.reset();
    m_timer_setup_solvercall.reset();
    m_timer_solve_assembly.reset();
    m_timer_solve_solvercall.reset();
}

}  // namespace chrono
