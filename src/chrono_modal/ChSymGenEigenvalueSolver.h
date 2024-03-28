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

#ifndef CHSYMGENEIGENVALUESOLVER_H
#define CHSYMGENEIGENVALUESOLVER_H

#include "chrono_modal/ChApiModal.h"
#include "chrono_modal/ChGeneralizedEigenvalueSolver.h"
#include "chrono/core/ChMatrix.h"
#include <Eigen/SparseCore>

#include <complex>

namespace chrono {

// Forward references
class ChDirectSolverLScomplex;
class ChDirectSolverLS;

namespace modal {

/// Base interface class for iterative eigenvalue solvers for generalized problem A*x = lambda*B*x
class ChApiModal ChSymGenEigenvalueSolver : public ChGeneralizedEigenvalueSolver<double> {
  public:
    ChSymGenEigenvalueSolver() {}
    virtual ~ChSymGenEigenvalueSolver(){};

    /// Solve the generalized eigenvalue problem A*x = lambda*B*x
    virtual int Solve(
        const ChSparseMatrix& A,
        const ChSparseMatrix& B,
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with n eigenvalues, will be resized
        int num_modes,
        ScalarType shift) const = 0;

    static void GetNaturalFrequencies(const ChVectorDynamic<ScalarType>& eigvals, ChVectorDynamic<double>& freq);

    static double GetNaturalFrequency(ScalarType eigval) { return sqrt(std::abs(eigval)) / CH_2PI; };

    static ScalarType GetOptimalShift(double freq) { return -std::pow(freq * CH_2PI, 2); }

};

/// Generalized eigenvalue solver implementing Krylov-Schur iterative method.
/// - generalized eigenvalue problem
/// - symmetric real sparse matrices
/// - shift-and-invert + real shift
class ChApiModal ChSymGenEigenvalueSolverKrylovSchur : public ChSymGenEigenvalueSolver {
  public:
    ChSymGenEigenvalueSolverKrylovSchur() {}
    virtual ~ChSymGenEigenvalueSolverKrylovSchur(){};

    /// Solve the eigenvalue problem
    virtual int Solve(
        const ChSparseMatrix& A,
        const ChSparseMatrix& B,
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with n eigenvalues, will be resized.
        int num_modes,
        ScalarType shift) const override;
};

/// Generalized eigenvalue solver implementing Lanczos iterative method.
/// - generalized eigenvalue problem
/// - symmetric real sparse matrices
/// - shift-and-invert + real shift
class ChApiModal ChSymGenEigenvalueSolverLanczos : public ChSymGenEigenvalueSolver {
  public:
    ChSymGenEigenvalueSolverLanczos() {}
    virtual ~ChSymGenEigenvalueSolverLanczos(){};

    /// Solve the eigenvalue problem
    virtual int Solve(
        const ChSparseMatrix& A,
        const ChSparseMatrix& B,
        ChMatrixDynamic<ScalarType>& eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
        ChVectorDynamic<ScalarType>& eigvals,   ///< output vector with n eigenvalues, will be resized.
        int num_modes,
        ScalarType shift) const override;
};

}  // end namespace modal

}  // end namespace chrono

#endif
