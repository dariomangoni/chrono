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

#ifndef CHSOLVEROSQP_H
#define CHSOLVEROSQP_H

#include "chrono_osqp/ChApiOsqp.h"
#include "chrono/solver/ChIterativeSolverVI.h"

#include <osqp.h>

namespace chrono {

/// @addtogroup osqp_module
/// @{

/** \class ChSolverOSQP
\brief Interface to the OSQP solver.


See ChIterativeSolverVI for more details.

<div class="ce-warning">
If appropriate and warranted by the problem setup, it is \e highly recommended to enable the sparsity pattern \e lock.
This can significantly improve performance for more complex problems (larger size and/or problems which include
constraints).
</div>

Minimal usage example, to be put anywhere in the code, before starting the main simulation loop:
\code{.cpp}
auto osqp_solver = chrono_types::make_shared<ChSolverOSQP>();
system.SetSolver(osqp_solver);
\endcode

See ChSystemDescriptor for more information about the problem formulation and the data structures passed to the solver.
*/
class ChApiOsqp ChSolverOSQP : public ChIterativeSolverVI {

    public:


    int runTest(){
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
            //data = (OSQPData *)c_malloc(sizeof(OSQPData));
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
        osqp_update_lin_cost(work, q_new);
        osqp_update_bounds(work, l_new, u_new);

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
    };

    // Inherited via ChIterativeSolverVI
    virtual double Solve(ChSystemDescriptor& sysd) override;
    virtual double GetError() const override;
};

/// @} osqp_module

}  // end namespace chrono

#endif
