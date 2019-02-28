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

#ifndef CHCVXOPTENGINE_H
#define CHCVXOPTENGINE_H

#include "chrono/core/ChCSMatrix.h"
#include "chrono_cvxopt/ChApiCvxopt.h"

#ifdef _DEBUG
#undef _DEBUG
#include <cvxopt.h>
#define _DEBUG
#else
#include <cvxopt.h>
#endif
#include "chrono_cvxopt/ChCSMatrixIndexType.h"

namespace chrono {
/// \class ChCvxoptConeQpEngine
/// Class that wraps the CVXOPT ConeQp solver.
/// It can solve VI and complementarity problems.

class ChApiCvxopt ChCvxoptConeQpEngine {
  public:
    explicit ChCvxoptConeQpEngine();
    virtual ~ChCvxoptConeQpEngine();

    void SetPMatrix(ChCSMatrixIndexType<int_t>& Pmat);

    void SetGMatrix(ChCSMatrixIndexType<int_t>& Gmat);

    void SetAMatrix(ChCSMatrixIndexType<int_t>& Amat);

    void SetqMatrix(ChMatrix<>& qmat);

    void SethMatrix(ChMatrix<>& hmat);

    void SetbMatrix(ChMatrix<>& bmat);

    /// Inform CVXOPT on the total number of unilateral linear constraints.
    void SetNumConstrUnilateralLinear(int_t constr_unilin_dim);

    /// Inform CVXOPT on the dimension of each cone constraint.
    /// \a unicon_dims must have an element for each cone constraint;
    /// \a unicon_dims[i] element reports the dimension of the i-th cone.
    void SetNumConstrUnilateralConic(const std::vector<int_t>& unicon_dims);

    void GetSolution(ChMatrix<>& x) const;

    void GetSlackVariable(ChMatrix<>& s) const;

    void GetLagrangianMultBilateral(ChMatrix<>& y) const;

    void GetLagrangianMultUnilateral(ChMatrix<>& z) const;

    // Returns true if the last ::Run returned an optimal solution;
    bool IsSolutionOptimal() const;

    bool SetMaxIters(int_t num = 100){ return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "maxiters", PyInt_FromSsize_t(num)); }
    bool SetAbsTol(double num = 1e-7) { return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "abstol", PyFloat_FromDouble(num)); }
    bool SetRelTol(double num = 1e-6) { return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "reltol", PyFloat_FromDouble(num)); }
    bool SetFeasTol(double num = 1e-7) { return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "feastol", PyFloat_FromDouble(num)); }
    bool SetRefinement(int_t num = 1) { return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "refinement", PyInt_FromSsize_t(num)); }

    int_t GetMaxIters() const { return PyInt_AsSsize_t(PyDict_GetItemString(PyObject_GetAttrString(solverModule, "options"), "maxiters")); }
    double GetAbsTol() const { return PyFloat_AsDouble(PyDict_GetItemString(PyObject_GetAttrString(solverModule, "options"), "abstol")); }
    double GetRelTol() const { return PyFloat_AsDouble(PyDict_GetItemString(PyObject_GetAttrString(solverModule, "options"), "reltol")); }
    double GetFeasTol() const { return PyFloat_AsDouble(PyDict_GetItemString(PyObject_GetAttrString(solverModule, "options"), "feastol")); }
    int_t GetRefinement() const { return PyInt_AsSsize_t(PyDict_GetItemString(PyObject_GetAttrString(solverModule, "options"), "refinement")); }

    double GetPrimalObjective() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "primal objective")); }
    double GetDualObjective() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "dual objective")); }
    double GetGap() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "gap")); }
    double GetRelativeGap() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "relative gap")); }
    double GetPrimalInfeasibility() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "primal infeasibility")); }
    double GetDualInfeasibility() const { return PyFloat_AsDouble(PyDict_GetItemString(sol, "dual infeasibility")); }

    void CheckMatrices() const;

    void Run();

    static void GetPyhtonError(bool print_fulltrace = false);

    bool SetShowProgress(bool sel) const { return !PyDict_SetItemString(PyObject_GetAttrString(solverModule, "options"), "show_progress", sel ? Py_True : Py_False); }

  private:
    PyObject* solverModule = nullptr;
    PyObject* programSolver = nullptr;
    PyObject* dims = PyDict_New();
    PyObject* sol = nullptr;

    spmatrix* P = nullptr;
    spmatrix* G = nullptr;
    spmatrix* A = nullptr;
    matrix* q = nullptr;
    matrix* h = nullptr;
    matrix* b = nullptr;

  protected:
    bool getDenseMatFromSolution(ChMatrix<>& z, std::string key) const;

};

}  // end of namespace chrono

#endif
