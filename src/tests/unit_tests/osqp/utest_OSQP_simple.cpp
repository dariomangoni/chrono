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
//
// ChronoOSQP unit test for OSQP
// =============================================================================

#include "gtest/gtest.h"
#include "chrono_osqp/ChOsqpEngine.h"

#include <osqp.h>

using namespace chrono;

const double ABS_ERR_D = 1e-8;

TEST(ChOsqpEngine, SetupSolve){

    ChVectorDynamicMap<double> x_Chrono_map(NULL, 0);
    ChVectorDynamicMap<double> y_Chrono_map(NULL, 0);

    double* x_Chrono;
    double* y_Chrono;

    double* x_Osqp;
    double* y_Osqp;

    const c_int n = 2;
    const c_int m = 3;

    {
        /// Data Loading
        ChSparseMatrixCSC P(n,n);
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletListP;
        tripletListP.reserve(3);

        tripletListP.push_back(T(0, 0, 4.0));
        tripletListP.push_back(T(0, 1, 1.0));
        tripletListP.push_back(T(1, 1, 2.0));

        P.setFromTriplets(tripletListP.begin(), tripletListP.end());

        ChSparseMatrixCSC A(m,n);
        std::vector<T> tripletListA;
        tripletListA.reserve(4);

        tripletListA.push_back(T(0, 0, 1.0));
        tripletListA.push_back(T(0, 1, 1.0));
        tripletListA.push_back(T(1, 0, 1.0));
        tripletListA.push_back(T(2, 1, 1.0));

        A.setFromTriplets(tripletListA.begin(), tripletListA.end());

        // Can be used ChVectorDynamic as well
        ChVectorN<double, static_cast<int>(n)> q = {1.0, 1.0};
        ChVectorN<double, m> l = {1.0, 0.0, 0.0};
        ChVectorN<double, m> u = {1.0, 0.7, 0.7};
    
        // OSQP Engine
        ChOsqpEngine osqp_engine;

        osqp_engine.SetProblem(P, A, q, l, u);
        osqp_engine.SetAlpha(1.0);
        ASSERT_FALSE(osqp_engine.Setup());
        ASSERT_FALSE(osqp_engine.Solve());

        // Extract solution
        osqp_engine.GetPrimalSolution(x_Chrono_map);
        osqp_engine.GetLagrangeMultipliers(y_Chrono_map);

        x_Chrono = osqp_engine.GetOsqpWorkspace()->solution->x;
        y_Chrono = osqp_engine.GetOsqpWorkspace()->solution->y;

        //std::cout << "PrimalSolution: ";
        //for(int x_sel=0; x_sel<P.rows(); ++x_sel){
        //    std::cout << x_Chrono(x_sel) << ", ";
        //}
        //std::cout << std::endl;

        //std::cout << "LagrangeMultipliers: ";
        //for(int y_sel=0; y_sel<A.rows(); ++y_sel){
        //    std::cout << y_Chrono(y_sel) << ", ";
        //}
        //std::cout << std::endl;

    }
    {

        // Load problem data
        c_float P_x[3] = {4.0, 1.0, 2.0, };
        c_int P_nnz = 3;
        c_int P_i[3] = {0, 0, 1, };
        c_int P_p[3] = {0, 1, 3, };
        c_float q[2] = {1.0, 1.0, };
        c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
        c_int A_nnz = 4;
        c_int A_i[4] = {0, 1, 0, 2, };
        c_int A_p[3] = {0, 2, 4, };
        c_float l[3] = {1.0, 0.0, 0.0, };
        c_float u[3] = {1.0, 0.7, 0.7, };


        // Exitflag
        c_int exitflag = 0;

        // Workspace structures
        OSQPWorkspace *work;
        OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

        // Populate data
        if (data) {
            data->n = n;
            data->m = m;
            data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
            data->q = q;
            data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
            data->l = l;
            data->u = u;
        }

        // Define solver settings as default
        if (settings) {
            osqp_set_default_settings(settings);
            settings->alpha = 1.0; // Change alpha parameter
        }

        // Setup workspace
        ASSERT_FALSE(osqp_setup(&work, data, settings));

        // Solve Problem
        ASSERT_FALSE(osqp_solve(work));

        // Cleanup
        if (data) {
            if (data->A) c_free(data->A);
            if (data->P) c_free(data->P);
            c_free(data);
        }
        if (settings) c_free(settings);

        x_Osqp = work->solution->x;
        y_Osqp = work->solution->y;

        //std::cout << "PrimalSolution: ";
        //for(int x_sel=0; x_sel<n; ++x_sel){
        //    std::cout << work->solution->x[x_sel] << ", ";
        //}
        //std::cout << std::endl;

        //std::cout << "LagrangeMultipliers: ";
        //for(int y_sel=0; y_sel<m; ++y_sel){
        //    std::cout << work->solution->y[y_sel] << ", ";
        //}
        //std::cout << std::endl;
    
    }

    for(int x_sel=0; x_sel<n; ++x_sel){
        ASSERT_NEAR(x_Chrono[x_sel], x_Chrono_map[x_sel], ABS_ERR_D);
    }

    for(int x_sel=0; x_sel<n; ++x_sel){
        ASSERT_NEAR(x_Osqp[x_sel], x_Chrono[x_sel], ABS_ERR_D);
    }

    for(int y_sel=0; y_sel<m; ++y_sel){
        ASSERT_NEAR(y_Chrono[y_sel], y_Chrono_map[y_sel], ABS_ERR_D);
    }

    for(int y_sel=0; y_sel<m; ++y_sel){
        ASSERT_NEAR(y_Osqp[y_sel], y_Chrono[y_sel], ABS_ERR_D);
    }



}


int testOsqpUpdateVectors(){
    // Load problem data
    c_float P_x[3] = {4.0, 1.0, 2.0, };
    c_int P_nnz = 3;
    c_int P_i[3] = {0, 0, 1, };
    c_int P_p[3] = {0, 1, 3, };
    c_float q[2] = {1.0, 1.0, };
    c_float q_new[2] = {2.0, 3.0, };
    c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
    c_int A_nnz = 4;
    c_int A_i[4] = {0, 1, 0, 2, };
    c_int A_p[3] = {0, 2, 4, };
    c_float l[3] = {1.0, 0.0, 0.0, };
    c_float l_new[3] = {2.0, -1.0, -1.0, };
    c_float u[3] = {1.0, 0.7, 0.7, };
    c_float u_new[3] = {2.0, 2.5, 2.5, };
    c_int n = 2;
    c_int m = 3;

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define solver settings as default
    if (settings) osqp_set_default_settings(settings);

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve problem
    osqp_solve(work);

    // Update problem
    osqp_update_lin_cost(work, q_new);
    osqp_update_bounds(work, l_new, u_new);

    // Solve updated problem
    osqp_solve(work);

    // Cleanup
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) c_free(settings);

    return exitflag;
}

int testOsqpUpdateMatrices(){
    // Load problem data
    c_float P_x[3] = {4.0, 1.0, 2.0, };
    c_float P_x_new[3] = {5.0, 1.5, 1.0, };
    c_int P_nnz = 3;
    c_int P_i[3] = {0, 0, 1, };
    c_int P_p[3] = {0, 1, 3, };
    c_float q[2] = {1.0, 1.0, };
    c_float q_new[2] = {2.0, 3.0, };
    c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
    c_float A_x_new[4] = {1.2, 1.5, 1.1, 0.8, };
    c_int A_nnz = 4;
    c_int A_i[4] = {0, 1, 0, 2, };
    c_int A_p[3] = {0, 2, 4, };
    c_float l[3] = {1.0, 0.0, 0.0, };
    c_float l_new[3] = {2.0, -1.0, -1.0, };
    c_float u[3] = {1.0, 0.7, 0.7, };
    c_float u_new[3] = {2.0, 2.5, 2.5, };
    c_int n = 2;
    c_int m = 3;

    // Exitflag
    c_int exitflag = 0;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
        data = (OSQPData *)c_malloc(sizeof(OSQPData));
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;
    }

    // Define Solver settings as default
    if (settings) osqp_set_default_settings(settings);

    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve problem
    osqp_solve(work);

    // Update problem
    // NB: Update only upper triangular part of P
    osqp_update_P(work, P_x_new, OSQP_NULL, 3);
    osqp_update_A(work, A_x_new, OSQP_NULL, 4);

    // Solve updated problem
    osqp_solve(work);

    // Cleanup
    if (data) {
        if (data->A) c_free(data->A);
        if (data->P) c_free(data->P);
        c_free(data);
    }
    if (settings) c_free(settings);

    return exitflag;
}
