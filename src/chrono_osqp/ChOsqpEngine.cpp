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


#include "chrono_osqp/ChOsqpEngine.h"

namespace chrono {

void ChOsqpEngine::SetProblem(const ChSparseMatrixCSC& P,
                              const ChSparseMatrixCSC& A,
                              ChVectorRef q,
                              ChVectorRef l,
                              ChVectorRef u) {


    data->m = A.rows();
    data->n = A.cols();
    
    data->P = csc_matrix(
        data->n, // First dimension
        data->n, // Second dimension
        P.nonZeros(), // Maximum number of nonzero elements
        const_cast<c_float*>(P.valuePtr()), // Vector of data
        const_cast<c_int*>(P.innerIndexPtr()), // Vector of row indices
        const_cast<c_int*>(P.outerIndexPtr())); // Vector of column pointers

    data->q = q.data();

    data->A = csc_matrix(
        data->m, // First dimension
        data->n, // Second dimension
        A.nonZeros(), // Maximum number of nonzero elements
        const_cast<c_float*>(A.valuePtr()), // Vector of data
        const_cast<c_int*>(A.innerIndexPtr()), // Vector of row indices
        const_cast<c_int*>(A.outerIndexPtr())); // Vector of column pointers

    data->l = l.data();
    data->u = u.data();

}

}  // namespace chrono
