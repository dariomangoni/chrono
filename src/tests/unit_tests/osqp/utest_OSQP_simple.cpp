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

#include "chrono_osqp/ChSolverOSQP.h"

int main(int argc, char **argv) {
    chrono::ChSolverOSQP mySolver;
    int res = mySolver.runTest();

    return res;
};