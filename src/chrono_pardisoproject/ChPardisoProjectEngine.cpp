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
// Authors: Dario Mangoni, Radu Serban
// =============================================================================


#include "chrono_pardisoproject/ChPardisoProjectEngine.h"

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);

namespace chrono {

ChPardisoProjectEngine::ChPardisoProjectEngine(pardisoproject_SYM symmetry) :symmetry(symmetry) {
    iparm[2]  = num_procs;
   
    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */
    solver = 0; /* use sparse direct solver */

    int mtype_int = symmetry;
    pardisoinit(pt,  &mtype_int, &solver, iparm, dparm, &error);

    printf("WARNING: PardisoProject suggests you to: \"Set environment OMP_NUM_THREADS to 1\"\n");

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
         //return 1;
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
}

ChPardisoProjectEngine::~ChPardisoProjectEngine() {
    PardisoProjectCall(pardisoproject_PHASE::END);
}

void ChPardisoProjectEngine::SetProblem(const ChSparseMatrix& Z, ChVectorRef rhs, ChVectorRef sol) {
    SetMatrix(Z);
    SetRhsVector(rhs);
    SetSolutionVector(sol);
}

void ChPardisoProjectEngine::SetMatrix(const ChSparseMatrix& Z) {
    this->SetMatrix(
        Z.rows(),
        const_cast<int*>(Z.outerIndexPtr()),
        const_cast<int*>(Z.innerIndexPtr()),
        const_cast<double*>(Z.valuePtr())
    );
}

void ChPardisoProjectEngine::SetMatrix(int n, int *ia, int *ja, double *a) {
    this->n = n;
    this->a = a;
    this->ia = ia;
    this->ja = ja;
}

void ChPardisoProjectEngine::SetMatrixSymmetry(pardisoproject_SYM symmetry) {
    this->symmetry = symmetry;
}

void ChPardisoProjectEngine::SetRhsVector(ChVectorRef b) {
    this->b = b.data();
}

void ChPardisoProjectEngine::SetRhsVector(double* b) {
    this->b = b;
}

void ChPardisoProjectEngine::SetSolutionVector(ChVectorRef x) {
    this->x = x.data();
}

void ChPardisoProjectEngine::SetSolutionVector(double* x) {
    this->x = x;
}


int ChPardisoProjectEngine::PardisoProjectCall(pardisoproject_PHASE phase) {
    int phase_int = phase;
    int mtype_int = symmetry;

    // TODO: WARNING: shifting indeces to 1-based whenever this phases are requested
    if (phase == pardisoproject_PHASE::ANALYZE || phase == pardisoproject_PHASE::ANALYZE_FACTORIZE || phase == pardisoproject_PHASE::COMPLETE)
        ShiftMatrixIndeces(+1);

    // TODO: assure to have non nullptr in appropriate arrays

    pardiso(pt, &maxfct, &mnum, &mtype_int, &phase_int, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    return error;
}

int ChPardisoProjectEngine::CheckMatrix(){
    /* -------------------------------------------------------------------- */
    /*  .. pardiso_chk_matrix(...)                                          */
    /*     Checks the consistency of the given matrix.                      */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    int mtype_int = symmetry;
    pardiso_chkmatrix(&mtype_int, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        return 1;
    }
    return error;
}

int ChPardisoProjectEngine::PrintStats() {
    /* -------------------------------------------------------------------- */
    /* .. pardiso_printstats(...)                                           */
    /*    prints information on the matrix to STDOUT.                       */
    /*    Use this functionality only for debugging purposes                */
    /* -------------------------------------------------------------------- */
    int mtype_int = symmetry;
    pardiso_printstats(&mtype_int, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        return 1;
    }
    return error;
}

void ChPardisoProjectEngine::shiftIndices(int n, int nonzeros, int* ia, int* ja, int value) {
    for (int i = 0; i < n + 1; i++) {
        ia[i] += value;
    }
    for (int i = 0; i < nonzeros; i++) {
        ja[i] += value;
    }
}

int ChPardisoProjectEngine::CheckVectors(){
    /* -------------------------------------------------------------------- */
    /* ..  pardiso_chkvec(...)                                              */
    /*     Checks the given vectors for infinite and NaN values             */
    /*     Input parameters (see PARDISO user manual for a description):    */
    /*     Use this functionality only for debugging purposes               */
    /* -------------------------------------------------------------------- */
    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        return 1;
    }
    return 0;
}

void ChPardisoProjectEngine::ShiftMatrixIndeces(int val){
    shiftIndices(n, ia[n], ia, ja, val); 
}

}  // namespace chrono
