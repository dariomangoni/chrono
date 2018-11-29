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
/*
A modified version of the example in the Python documentation:
http://docs.python.org/ext/pure-embedding.html

Solves the LP:

minimize 0.5*xT*P*x + qT*x
subject  G*x + s = h
A*x = b
s >= 0
where >= is intended as general inequality.

P=
2,14	-0,47	-2,33
-0,47	 4,87	-0,96
-2,33	-0,96	 5,89

G=
-1	0	0
0	-1	0
0	0	-1
0	0	0
1	0	0
0	1	0
0	0	1

q=
-0,97
-2,73
0,33

h=
0
0
0
1
0
0
0


Solution i.e. x is
[ 7.26e-01]
[ 6.18e-01]
[ 3.03e-01]

and prints the solution from C.


On Ubuntu Linux compile with:

gcc -o embed_cvxopt embed_cvxopt.c \
-I${CVXOPT_SRC}/C \
-I/usr/include/python2.4 -lpython2.4
*/


#include "chrono_cvxopt/ChCSMatrixIndexType.h"
#include "core/ChCSMatrix.h"
#include "core/ChMatrixNM.h"
#include "chrono_cvxopt/ChCvxoptConeQpEngine.h"
#include "core/ChMatrixDynamic.h"
#ifdef _DEBUG
#undef _DEBUG
#include <cvxopt.h>
#define _DEBUG
#else
#include <cvxopt.h>
#endif


int main(int argc, char* argv[]) {
    const int m = 7;
    const int n = 3;
    const int nnzP = 9;
    const int nnzG = 6;

    chrono::ChCSMatrixIndexType<int_t> chronoP(n, n, false, 9);
    chrono::ChCSMatrixIndexType<int_t> chronoG(m, n, false, 6);
    chrono::ChMatrixNM<double, n, 1> chrono_q;
    chrono::ChMatrixNM<double, m, 1> chrono_h;

    chronoP.SetElement(0, 0, +2.14);
    chronoP.SetElement(1, 0, -0.47);
    chronoP.SetElement(2, 0, -2.33);
    chronoP.SetElement(0, 1, -0.47);
    chronoP.SetElement(1, 1, +4.87);
    chronoP.SetElement(2, 1, -0.96);
    chronoP.SetElement(0, 2, -2.33);
    chronoP.SetElement(1, 2, -0.96);
    chronoP.SetElement(2, 2, +5.89);
    chronoP.Compress();

    chronoG.SetElement(0, 0, -1.0);
    chronoG.SetElement(4, 0, +1.0);
    chronoG.SetElement(1, 1, -1.0);
    chronoG.SetElement(5, 1, +1.0);
    chronoG.SetElement(2, 2, -1.0);
    chronoG.SetElement(6, 2, +1.0);
    chronoG.Compress();

    chrono_q.SetElement(0, 0, -0.97);
    chrono_q.SetElement(1, 0, -2.73);
    chrono_q.SetElement(2, 0, +0.33);

    chrono_h.SetElement(3, 0, 1.0);

    // Setup the CVXOPT engine to solve the problem
    chrono::ChCvxoptConeQpEngine engine;
    chrono::ChMatrixDynamic<double> sol, lagmul;

    engine.SetPMatrix(chronoP);
    engine.SetGMatrix(chronoG);
    engine.SetqMatrix(chrono_q);
    engine.SethMatrix(chrono_h);
    engine.SetNumConstrUnilateralLinear(3);
    engine.SetNumConstrUnilateralConic(std::vector<int_t>(1, 4));

    engine.CheckMatrices();

    engine.Run();

    engine.GetSolution(sol);
    engine.GetLagrangianMultUnilateral(lagmul);

    chrono::ChMatrixDynamic<double> expected_sol;
    expected_sol.Resize(3, 1);
    expected_sol.SetElement(0, 0, 0.72558318685981849);
    expected_sol.SetElement(1, 0, 0.61806264311119241);
    expected_sol.SetElement(2, 0, 0.30253527966423421);

    chrono::GetLog() << "Solution is:\n" << sol << "\n";
    chrono::GetLog() << "Lagrangian Multipliers is:\n" << lagmul << "\n";

    auto solution_reached = (sol - expected_sol).NormTwo()/sol.GetRows() < 1e-6;
    auto solution_optimal = engine.IsSolutionOptimal();

    return !(solution_reached && solution_optimal);

}