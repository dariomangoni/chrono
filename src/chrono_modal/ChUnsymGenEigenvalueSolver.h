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

#ifndef CHUNSYMGENEIGENVALUESOLVER_H
#define CHUNSYMGENEIGENVALUESOLVER_H

#include "chrono_modal/ChApiModal.h"
#include "chrono/solver/ChDirectSolverLScomplex.h"
#include "chrono_modal/ChGeneralizedEigenvalueSolver.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChTimer.h"
#include "chrono/physics/ChAssembly.h"

#include <complex>

namespace chrono {

// Forward references
class ChDirectSolverLScomplex;
class ChDirectSolverLS;

namespace modal {

/// Base interface class for generalized eigenvalue problems with real generic (i.e. potentially unsymmetric) matrices.
/// It is currently assumed that all the derived classes implement shift-and-invert strategy with a complex shift.
class ChApiModal ChUnsymGenEigenvalueSolver : public ChGeneralizedEigenvalueSolver<std::complex<double>> {
  public:
    ChUnsymGenEigenvalueSolver() {}
    virtual ~ChUnsymGenEigenvalueSolver(){};

    virtual int Solve(
        const ChSparseMatrix& A,                ///< input A matrix
        const ChSparseMatrix& B,                ///< input B matrix
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with eigenvalues (real part
                                                ///< not zero if some damping), will be resized
        int num_modes,    ///< optional: number of desired lower eigenvalues. If zero return all eigenvalues
        ScalarType sigma  ///< optional: shift for the shift-and-invert strategy
    ) const = 0;

    static void GetNaturalFrequencies(const ChVectorDynamic<ScalarType>& eigvals,
                                      ChVectorDynamic<double>& natural_freq);

    static void GetDampedFrequencies(const ChVectorDynamic<ScalarType>& eigvals, ChVectorDynamic<double>& damped_freq);

    static void GetDampingRatios(const ChVectorDynamic<ScalarType>& eigvals, ChVectorDynamic<double>& damp_ratios);

    static double GetNaturalFrequency(ScalarType eigval) { return std::abs(eigval) / CH_2PI; };

    static ScalarType GetOptimalShift(double frequency) { return frequency * CH_2PI; };
};

/// Solves the generalized eigenvalue problem with the Krylov-Schur iterative method.
/// Features:
/// - generalized eigenvalue problem
/// - generic sparse matrices
/// - shift-and-invert + complex shift
class ChApiModal ChUnsymGenEigenvalueSolverKrylovSchur : public ChUnsymGenEigenvalueSolver {
  public:
    /// Default: uses Eigen::SparseQR as factorization for the shift&invert,
    /// otherwise pass a custom complex sparse solver for faster factorization (ex. ChSolverComplexPardisoMKL)
    ChUnsymGenEigenvalueSolverKrylovSchur(
        std::shared_ptr<ChDirectSolverLScomplex> linear_solver = chrono_types::make_shared<ChSolverSparseComplexLU>());

    virtual ~ChUnsymGenEigenvalueSolverKrylovSchur(){};

    /// Solve the quadratic eigenvalue problem (lambda^2*M + lambda*R + K)*x = 0 s.t. Cq*x = 0
    /// If n_modes=0, return all eigenvalues, otherwise only the first lower n_modes.
    virtual int Solve(
        const ChSparseMatrix& A,                ///< input A matrix
        const ChSparseMatrix& B,                ///< input B matrix
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with eigenvalues (real part
                                                ///< not zero if some damping), will be resized
        int num_modes,    ///< optional: number of desired lower eigenvalues. If zero return all eigenvalues
        ScalarType sigma  ///< optional: shift for the shift-and-invert strategy
    ) const override;

  protected:
    std::shared_ptr<ChDirectSolverLScomplex> m_linear_solver;
};

}  // end namespace modal

}  // end namespace chrono

#endif
