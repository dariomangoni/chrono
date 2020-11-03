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

#ifndef CHPARDISOPROJECTENGINE_H
#define CHPARDISOPROJECTENGINE_H

#include "chrono/core/ChMatrix.h"
#include "chrono_pardisoproject/ChApiPardisoProject.h"


namespace chrono {

/// @addtogroup pardisoproject_module
/// @{

/// Wrapper class for the PardisoProject direct linear solver.
/// This solver is not appropriate for VI and complementarity problems.
class ChApiPardisoProject ChPardisoProjectEngine {
  public:
    enum pardisoproject_SYM {
        STRUCTURAL_SYMMETRIC = 1, ///< real and structurally symmetric, supernode pivoting
        SYMMETRIC_POSDEF = 2, ///< real and symmetric positive definite
        SYMMETRIC_GENERAL = -2, ///< real and symmetric indefinite, diagonal or Bunch-Kaufman pivoting
        UNSYMMETRIC = 11 ///< real and nonsymmetric, complete supernode pivoting
    };

    enum pardisoproject_PHASE {
        END = -1,
        ANALYZE = 11,
        ANALYZE_FACTORIZE = 12,
        FACTORIZE = 22,
        SOLVE = 33,
        FACTORIZE_SOLVE = 23,
        COMPLETE = 13,
        SELECTED_INVERSION = -22
    };

    ChPardisoProjectEngine(pardisoproject_SYM symmetry);
    ~ChPardisoProjectEngine();

    /// Set the problem matrix and the right-hand side.
    void SetProblem(const ChSparseMatrix& Z, ChVectorRef rhs, ChVectorRef sol);

    /// Set the problem matrix.
    void SetMatrix(const ChSparseMatrix& Z);
    void SetMatrix(int n, int *ia, int *ja, double *a);

    /// Informs PardisoProject of the matrix symmetry type.
    void SetMatrixSymmetry(pardisoproject_SYM symmetry);

    /// Set the right-hand side vector.
    /// Note that it is the caller's responsibility to ensure that the size is appropriate.
    void SetRhsVector(ChVectorRef b);
    void SetRhsVector(double* b);

    /// Set the solution vector.
    /// Note that it is the caller's responsibility to ensure that the size is appropriate.
    void SetSolutionVector(ChVectorRef x);
    void SetSolutionVector(double* x);

    /// Submit job to PardisoProject.
    int PardisoProjectCall(pardisoproject_PHASE job_call);

    int CheckMatrix();

    int CheckVectors();

    int PrintStats();

    int GetIPARM(int id) const { return iparm[id]; };
    void SetIPARM(int id, int val){ iparm[id] = val; };

    int GetLastError() {return error;};

    void ShiftMatrixIndeces(int val);


  private:
    void* m_engine;  ///< underlying PardisoProject interface

    /* RHS and solution vectors. */
    double   *b, *x;
    int      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      solver;
    int      maxfct, mnum, phase, error, msglvl;

    /* Number of processors. */
    int      num_procs = 1;

    /* Auxiliary variables. */
    char    *var;
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */
    
    // Matrix variables
    int    n = 0;
    int    *ia;
    int    *ja;
    double  *a;

    pardisoproject_SYM symmetry;

    void shiftIndices(int n, int nonzeros, int* ia, int* ja, int value);


};

/// @} pardisoproject_module

}  // end of namespace chrono

#endif
