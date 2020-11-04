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

#include "chrono_pardisoproject/ChSolverPardisoProject.h"
#include "chrono/parallel/ChOpenMP.h"


namespace chrono {
ChSolverPardisoProject::ChSolverPardisoProject(int num_threads, ChPardisoProjectEngine::pardisoproject_SYM symmetry)
    : m_engine(symmetry) {
    //int nthreads = (num_threads <= 0) ? ChOMP::GetNumProcs() : num_threads;
    //ChOMP::SetNumThreads(nthreads);
    m_engine.SetIPARM(5,0);
}


bool ChSolverPardisoProject::FactorizeMatrix() {
    m_engine.SetMatrix(m_mat);

    m_engine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::ANALYZE);
    m_engine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::FACTORIZE);

    
    if (verbose){
        if (m_engine.GetLastError() != 0) {
        printf("\nERROR during symbolic factorization: %d", m_engine.GetLastError());
            exit(1);
        }
        printf("\nReordering completed ... ");
        printf("\nNumber of nonzeros in factors  = %d", m_engine.GetIPARM(17));
        printf("\nNumber of factorization MFLOPS = %d", m_engine.GetIPARM(18));
    }


    return true;
}

bool ChSolverPardisoProject::SolveSystem() {
    m_engine.SetRhsVector(m_rhs);
    m_engine.SetSolutionVector(m_sol);
    m_engine.PardisoProjectCall(ChPardisoProjectEngine::pardisoproject_PHASE::SOLVE);

   
    if (m_engine.GetLastError() != 0) {
        printf("\nERROR during solution: %d", m_engine.GetLastError());
        exit(3);
    }
    return true;
}

void ChSolverPardisoProject::PrintErrorMessage() {
    printf("\nERROR: %d", m_engine.GetLastError());
}



}  // end of namespace chrono
