/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on for unsymmetric linear systems                               */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */

// C++ compatible

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono_pardisoproject\ChPardisoProjectEngine.h>

using namespace std;
using namespace chrono;

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);


int main( void ) 
{
    ChPardisoProjectEngine ppengine(ChPardisoProjectEngine::pardisoproject_SYM::UNSYMMETRIC);
    /* Matrix data. */
    int    n = 8;
    int    ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
    int    ja[20] = { 0,    2,       5, 6, 
                         1, 2,    4,
                            2,             7,
                               3,       6,
                         1,
                            2,       5,    7,
                         1,             6,
                            2,          6, 7 };
    double  a[20] = { 7.0,      1.0,           2.0, 7.0, 
                          -4.0, 8.0,      2.0,
                                1.0,                     5.0,
                                     7.0,           9.0,
                          -4.0,
                                7.0,           3.0,      8.0,
                           1.0,                    11.0,
                               -3.0,                2.0, 5.0 };

    int      nnz = ia[n];
    int      mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */
    double   b[8], x[8];

    ppengine.SetMatrix(n, ia, ja, a);
    ppengine.SetSolutionVector(x);
    ppengine.SetRhsVector(b);

//
///* -------------------------------------------------------------------- */
///*  .. pardiso_chk_matrix(...)                                          */
///*     Checks the consistency of the given matrix.                      */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//    
//    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//    if (error != 0) {
//        printf("\nERROR in consistency of matrix: %d", error);
//        exit(1);
//    }
//
///* -------------------------------------------------------------------- */
///* ..  pardiso_chkvec(...)                                              */
///*     Checks the given vectors for infinite and NaN values             */
///*     Input parameters (see PARDISO user manual for a description):    */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//
//    pardiso_chkvec (&n, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR  in right hand side: %d", error);
//        exit(1);
//    }
//
///* -------------------------------------------------------------------- */
///* .. pardiso_printstats(...)                                           */
///*    prints information on the matrix to STDOUT.                       */
///*    Use this functionality only for debugging purposes                */
///* -------------------------------------------------------------------- */
//
//    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
//    if (error != 0) {
//        printf("\nERROR right hand side: %d", error);
//        exit(1);
//    }
// 
/* -------------------------------------------------------------------- */    
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/*     ppengine.PardisoProjectCall
 */ 
    ppengine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::ANALYZE);
    
    if (ppengine.GetLastError() != 0) {
        printf("\nERROR during symbolic factorization: %d", ppengine.GetLastError());
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", ppengine.GetIPARM(17));
    printf("\nNumber of factorization MFLOPS = %d", ppengine.GetIPARM(18));
   
/* -------------------------------------------------------------------- */    
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    
    ppengine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::FACTORIZE);

    if (ppengine.GetLastError() != 0) {
        printf("\nERROR during numerical factorization: %d", ppengine.GetLastError());
        exit(2);
    }
    printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    ppengine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::SOLVE);

   
    if (ppengine.GetLastError() != 0) {
        printf("\nERROR during solution: %d", ppengine.GetLastError());
        exit(3);
    }

    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (int i = 0; i < n; i++) {
        printf("\n x [%d] = % f", i, x[i] );
    }
    printf ("\n");

///* -------------------------------------------------------------------- */
///* ..  Back substitution with tranposed matrix A^t x=b                   */
///* -------------------------------------------------------------------- */
//    phase = 33;
//
//    iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
//    iparm[11] = 1;       /* Solving with transpose matrix. */
//
//    /* Set right hand side to one. */
//    for (i = 0; i < n; i++) {
//        b[i] = 1;
//    }
//  
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, b, x, &error,  dparm);
//  
//    if (error != 0) {
//        printf("\nERROR during solution: %d", error);
//        exit(3);
//    }
//
//    printf("\nSolve completed ... ");
//    printf("\nThe solution of the system is: ");
//    for (i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, x[i] );
//    }
//    printf ("\n");

///* -------------------------------------------------------------------- */    
///* ..  Convert matrix back to 0-based C-notation.                       */
///* -------------------------------------------------------------------- */ 
//    for (i = 0; i < n+1; i++) {
//        ia[i] -= 1;
//    }
//    for (i = 0; i < nnz; i++) {
//        ja[i] -= 1;
//    }
    ppengine.ShiftMatrixIndeces(-1);

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */ 

    return 0;
} 