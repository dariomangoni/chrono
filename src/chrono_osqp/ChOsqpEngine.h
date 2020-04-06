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

#ifndef CHOSQPENGINE_H
#define CHOSQPENGINE_H

#include "chrono/core/ChMatrix.h"
#include "chrono_osqp/ChApiOsqp.h"

#include <osqp.h>



namespace chrono {


using ChSparseMatrixCSC = Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int>;

template <typename T>
using ChVectorDynamicMap = Eigen::Map<ChVectorDynamic<T>>;

/// @addtogroup osqp_module
/// @{

/// Wrapper class for the OSQP solver.
/// This solver is appropriate for VI and complementarity problems.
class ChApiOsqp ChOsqpEngine {
  public:

    ChOsqpEngine(){
        settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        data = (OSQPData *)c_malloc(sizeof(OSQPData));
        osqp_set_default_settings(settings);
    };
    ~ChOsqpEngine(){
        c_free(data); // do not free internal P and A arrays!
    };

    /// Set the problem matrix and the right-hand side.
    void SetProblem(const ChSparseMatrixCSC& P, const ChSparseMatrixCSC& A, ChVectorRef q, ChVectorRef l, ChVectorRef u);

    c_int Setup() {
        return osqp_setup(&work, data, settings);
    };

    c_int Solve() {
        return osqp_solve(work);
    };

    c_int UpdateVectors(ChVectorRef q_new, ChVectorRef l_new, ChVectorRef u_new) {
        c_int up_lin_cost_exitflag = osqp_update_lin_cost(work, q_new.data());
        c_int up_bound_exitflag = osqp_update_bounds(work, l_new.data(), u_new.data());
        return up_lin_cost_exitflag==0 && up_bound_exitflag==0 ? 0 : 1;
    };

    c_int UpdateMatrices(ChVectorRef P_x_new, ChVectorRef A_x_new) {
        c_int up_P_exitflag = osqp_update_P(work, P_x_new.data(), OSQP_NULL, 3);
        c_int up_A_exitflag = osqp_update_A(work, A_x_new.data(), OSQP_NULL, 4);
        return up_P_exitflag==0 && up_A_exitflag==0 ? 0 : 1;
    };

    void GetPrimalSolution(ChVectorDynamicMap<double>& x){
        ::new (&x) ChVectorDynamicMap<double>(work->solution->x, data->n);
    }

    void GetLagrangeMultipliers(ChVectorDynamicMap<double>& y){
        ::new (&y) ChVectorDynamicMap<double>(work->solution->y, data->m);
    }

    void SetAlpha(double alpha){
        settings->alpha = alpha;
    }


    OSQPWorkspace* GetOsqpWorkspace() { return work; }
    OSQPSettings* GetOsqpSettings() { return settings; }
    OSQPData* GetOswpData() {return data; }

  private:

    OSQPWorkspace* work = nullptr;
    OSQPSettings* settings;
    OSQPData* data;


};

/// @} osqp_module

}  // end of namespace chrono

#endif
