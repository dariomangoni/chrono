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

template <typename index_type = int_t>
class ChCSMatrix_int_t {
  public:
    explicit ChCSMatrix_int_t(chrono::ChCSMatrix& mat) : mat(mat) { CloneFromChCSMatrix(mat); }

    // void* values = nullptr;       /* value list */
    // index_type* colptr = nullptr; /* column pointer list */
    // index_type* rowind = nullptr; /* row index list */
    index_type nrows, ncols; /* number of rows and columns */
    int id = DOUBLE;         /* DOUBLE, COMPLEX */

    chrono::ChCSMatrix& mat;
    std::vector<index_type> leadIndex;
    std::vector<index_type> trailIndex;

    void CloneFromChCSMatrix(chrono::ChCSMatrix& mat) {
        nrows = mat.GetNumRows();
        ncols = mat.GetNumColumns();

        mat.Compress();

        leadIndex.resize(ncols + 1);
        trailIndex.resize(mat.GetNNZ());

        auto mat_trailIndx = mat.GetCS_TrailingIndexArray();
        auto mat_leadIndx = mat.GetCS_LeadingIndexArray();
        auto mat_values = mat.GetCS_ValueArray();

        for (auto col_sel = 0; col_sel < ncols + 1; ++col_sel) {
            leadIndex[col_sel] = mat_leadIndx[col_sel];
        }

        for (auto nnz_sel = 0; nnz_sel < mat.GetNNZ(); ++nnz_sel) {
            trailIndex[nnz_sel] = mat_trailIndx[nnz_sel];
        }
    }

    index_type* GetCS_TrailingIndexArray() { return trailIndex.data(); }
    index_type* GetCS_LeadingIndexArray() { return leadIndex.data(); }
    double* GetCS_ValueArray() { return values.data(); }
};

int demo_withMinimalMatrix() {
    const int m = 7;
    const int n = 3;
    const int nnzP = 9;
    const int nnzG = 6;

    chrono::ChCSMatrix chronoP(n, n, false, 9), chronoG(m, n, false, 6);
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

    /////////////////
    Py_Initialize();

    if (import_cvxopt() < 0) {
        fprintf(stderr, "error importing cvxopt");
        getchar();
        return 1;
    }

    /* import cvxopt.solvers */
    PyObject* solverModule = PyImport_ImportModule("cvxopt.solvers");
    if (!solverModule) {
        fprintf(stderr, "error importing cvxopt.solvers");
        getchar();
        return 1;
    }

    /* get reference to solvers.solvelp */
    PyObject* programSolver = PyObject_GetAttrString(solverModule, "coneqp");
    if (!programSolver) {
        fprintf(stderr, "error referencing cvxopt.solvers.coneqp");
        Py_DECREF(solverModule);
        getchar();
        return 1;
    }

    PyObject* P = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    PyObject* G = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    // PyObject* P = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    // PyObject* G = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    PyObject* q = (PyObject*)Matrix_New(0, 1, DOUBLE);
    PyObject* h = (PyObject*)Matrix_New(0, 1, DOUBLE);

    PyObject* pArgs = PyTuple_New(5);
    PyObject* dims = PyDict_New();
    if (!P || !q || !G || !h || !pArgs) {
        fprintf(stderr, "error creating matrices");
        Py_DECREF(solverModule);
        Py_DECREF(programSolver);
        Py_XDECREF(P);
        Py_XDECREF(q);
        Py_XDECREF(G);
        Py_XDECREF(h);
        Py_XDECREF(pArgs);
        getchar();
        return 1;
    }

    std::cout << std::endl << "Matrix BEFORE re-assignment" << std::endl;
    std::cout << "P: " << (SpMatrix_Check(P) ? "OK" : "ERROR") << std::endl;
    std::cout << "G: " << (SpMatrix_Check(G) ? "OK" : "ERROR") << std::endl;
    std::cout << "q: " << (Matrix_Check(q) ? "OK" : "ERROR") << std::endl;
    std::cout << "h: " << (Matrix_Check(h) ? "OK" : "ERROR") << std::endl;

    ChCSMatrix_int_t<int_t> P_int_t(chronoP), G_int_t(chronoG);

    Py_INCREF(P);
    auto Pccs_ptr = reinterpret_cast<spmatrix*>(P)->obj;
    Pccs_ptr->colptr = P_int_t.GetCS_LeadingIndexArray();
    Pccs_ptr->rowind = P_int_t.GetCS_TrailingIndexArray();
    Pccs_ptr->values = chronoP.GetCS_ValueArray();
    Pccs_ptr->nrows = P_int_t.nrows;
    Pccs_ptr->ncols = P_int_t.ncols;

    Py_INCREF(G);
    auto Gccs_ptr = reinterpret_cast<spmatrix*>(G)->obj;
    Gccs_ptr->colptr = G_int_t.GetCS_LeadingIndexArray();
    Gccs_ptr->rowind = G_int_t.GetCS_TrailingIndexArray();
    Gccs_ptr->values = chronoG.GetCS_ValueArray();
    Gccs_ptr->nrows = G_int_t.nrows;
    Gccs_ptr->ncols = G_int_t.ncols;

    Py_INCREF(q);
    auto q_ccs_ptr = reinterpret_cast<matrix*>(q);
    q_ccs_ptr->nrows = chrono_q.GetRows();
    q_ccs_ptr->ncols = chrono_q.GetColumns();
    q_ccs_ptr->buffer = chrono_q.GetAddress();

    Py_INCREF(h);
    auto h_ccs_ptr = reinterpret_cast<matrix*>(h);
    h_ccs_ptr->nrows = chrono_h.GetRows();
    h_ccs_ptr->ncols = chrono_h.GetColumns();
    h_ccs_ptr->buffer = chrono_h.GetAddress();

    std::cout << std::endl << "Matrix AFTER re-assignment" << std::endl;
    std::cout << "P: " << (SpMatrix_Check(P) ? "OK" : "ERROR") << std::endl;
    std::cout << "G: " << (SpMatrix_Check(G) ? "OK" : "ERROR") << std::endl;
    std::cout << "q: " << (Matrix_Check(q) ? "OK" : "ERROR") << std::endl;
    std::cout << "h: " << (Matrix_Check(h) ? "OK" : "ERROR") << std::endl;

    PyDict_SetItemString(dims, "l", PyInt_FromLong(n));
    PyObject* a = PyList_New(1);
    PyList_SetItem(a, 0, PyInt_FromLong(n + 1));
    PyDict_SetItemString(dims, "q", a);
    PyDict_SetItemString(dims, "s", PyList_New(0));

    /* pack matrices into an argument tuple - references are stolen*/
    PyTuple_SetItem(pArgs, 0, P);
    PyTuple_SetItem(pArgs, 1, q);
    PyTuple_SetItem(pArgs, 2, G);
    PyTuple_SetItem(pArgs, 3, h);
    PyTuple_SetItem(pArgs, 4, dims);

    PyObject* sol = PyObject_CallObject(programSolver, pArgs);
    if (!sol) {
        PyErr_Print();
        Py_DECREF(solverModule);
        Py_DECREF(programSolver);
        Py_DECREF(pArgs);
        getchar();
        return 1;
    }

    PyObject* x = PyDict_GetItemString(sol, "x");

    std::cout << "x: " << std::endl;
    for (auto row = 0; row < n; ++row)
        std::cout << MAT_BUFD(x)[row] << std::endl;

    Py_DECREF(solverModule);
    Py_DECREF(programSolver);
    Py_DECREF(pArgs);
    Py_DECREF(sol);

    Py_Finalize();

    getchar();
    return 0;
}

int demo_withChCSMatrixIndexType()
{
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
    chronoP.VerifyMatrix();

    chronoG.SetElement(0, 0, -1.0);
    chronoG.SetElement(4, 0, +1.0);
    chronoG.SetElement(1, 1, -1.0);
    chronoG.SetElement(5, 1, +1.0);
    chronoG.SetElement(2, 2, -1.0);
    chronoG.SetElement(6, 2, +1.0);
    chronoG.Compress();
    chronoP.VerifyMatrix();

    chrono_q.SetElement(0, 0, -0.97);
    chrono_q.SetElement(1, 0, -2.73);
    chrono_q.SetElement(2, 0, +0.33);

    chrono_h.SetElement(3, 0, 1.0);

    /////////////////
    Py_Initialize();

    if (import_cvxopt() < 0) {
        fprintf(stderr, "error importing cvxopt");
        getchar();
        return 1;
    }

    /* import cvxopt.solvers */
    PyObject* solverModule = PyImport_ImportModule("cvxopt.solvers");
    if (!solverModule) {
        fprintf(stderr, "error importing cvxopt.solvers");
        getchar();
        return 1;
    }

    /* get reference to solvers.solvelp */
    PyObject* programSolver = PyObject_GetAttrString(solverModule, "coneqp");
    if (!programSolver) {
        fprintf(stderr, "error referencing cvxopt.solvers.coneqp");
        Py_DECREF(solverModule);
        getchar();
        return 1;
    }

    PyObject* P = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    PyObject* G = (PyObject*)SpMatrix_New(0, 0, 0, DOUBLE);
    // PyObject* P = (PyObject*)SpMatrix_New(n, n, 0, DOUBLE);
    // PyObject* G = (PyObject*)SpMatrix_New(m, n, 0, DOUBLE);
    PyObject* q = (PyObject*)Matrix_New(0, 1, DOUBLE);
    PyObject* h = (PyObject*)Matrix_New(0, 1, DOUBLE);

    PyObject* pArgs = PyTuple_New(5);
    PyObject* dims = PyDict_New();
    if (!P || !q || !G || !h || !pArgs) {
        fprintf(stderr, "error creating matrices");
        Py_DECREF(solverModule);
        Py_DECREF(programSolver);
        Py_XDECREF(P);
        Py_XDECREF(q);
        Py_XDECREF(G);
        Py_XDECREF(h);
        Py_XDECREF(pArgs);
        getchar();
        return 1;
    }

    std::cout << std::endl << "Matrix BEFORE re-assignment" << std::endl;
    std::cout << "P: " << (SpMatrix_Check(P) ? "OK" : "ERROR") << std::endl;
    std::cout << "G: " << (SpMatrix_Check(G) ? "OK" : "ERROR") << std::endl;
    std::cout << "q: " << (Matrix_Check(q) ? "OK" : "ERROR") << std::endl;
    std::cout << "h: " << (Matrix_Check(h) ? "OK" : "ERROR") << std::endl;

    Py_INCREF(P);
    auto Pccs_ptr = reinterpret_cast<spmatrix*>(P)->obj;
    Pccs_ptr->colptr = chronoP.GetCS_LeadingIndexArray();
    Pccs_ptr->rowind = chronoP.GetCS_TrailingIndexArray();
    Pccs_ptr->values = reinterpret_cast<void*>(chronoP.GetCS_ValueArray());
    Pccs_ptr->nrows = chronoP.GetNumRows();
    Pccs_ptr->ncols = chronoP.GetNumColumns();

    Py_INCREF(G);
    auto Gccs_ptr = reinterpret_cast<spmatrix*>(G)->obj;
    Gccs_ptr->colptr = chronoG.GetCS_LeadingIndexArray();
    Gccs_ptr->rowind = chronoG.GetCS_TrailingIndexArray();
    Gccs_ptr->values = chronoG.GetCS_ValueArray();
    Gccs_ptr->nrows = chronoG.GetNumRows();
    Gccs_ptr->ncols = chronoG.GetNumColumns();

    Py_INCREF(q);
    auto q_ccs_ptr = reinterpret_cast<matrix*>(q);
    q_ccs_ptr->nrows = chrono_q.GetRows();
    q_ccs_ptr->ncols = chrono_q.GetColumns();
    q_ccs_ptr->buffer = chrono_q.GetAddress();

    Py_INCREF(h);
    auto h_ccs_ptr = reinterpret_cast<matrix*>(h);
    h_ccs_ptr->nrows = chrono_h.GetRows();
    h_ccs_ptr->ncols = chrono_h.GetColumns();
    h_ccs_ptr->buffer = chrono_h.GetAddress();

    std::cout << std::endl << "Matrix AFTER re-assignment" << std::endl;
    std::cout << "P: " << (SpMatrix_Check(P) ? "OK" : "ERROR") << std::endl;
    std::cout << "G: " << (SpMatrix_Check(G) ? "OK" : "ERROR") << std::endl;
    std::cout << "q: " << (Matrix_Check(q) ? "OK" : "ERROR") << std::endl;
    std::cout << "h: " << (Matrix_Check(h) ? "OK" : "ERROR") << std::endl;

    PyDict_SetItemString(dims, "l", PyInt_FromLong(n));
    PyObject* a = PyList_New(1);
    PyList_SetItem(a, 0, PyInt_FromLong(n + 1));
    PyDict_SetItemString(dims, "q", a);
    PyDict_SetItemString(dims, "s", PyList_New(0));

    /* pack matrices into an argument tuple - references are stolen*/
    PyTuple_SetItem(pArgs, 0, P);
    PyTuple_SetItem(pArgs, 1, q);
    PyTuple_SetItem(pArgs, 2, G);
    PyTuple_SetItem(pArgs, 3, h);
    PyTuple_SetItem(pArgs, 4, dims);

    PyObject* sol = PyObject_CallObject(programSolver, pArgs);
    if (!sol) {
        PyErr_Print();
        Py_DECREF(solverModule);
        Py_DECREF(programSolver);
        Py_DECREF(pArgs);
        getchar();
        return 1;
    }

    PyObject* x = PyDict_GetItemString(sol, "x");

    std::cout << "x: " << std::endl;
    for (auto row = 0; row < n; ++row)
        std::cout << MAT_BUFD(x)[row] << std::endl;

    Py_DECREF(solverModule);
    Py_DECREF(programSolver);
    Py_DECREF(pArgs);
    Py_DECREF(sol);

    Py_Finalize();

    getchar();
    return 0;
}

int demo_withEngine()
{
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

    /////////////////
    chrono::ChCvxoptConeQpEngine engine;
    chrono::ChMatrixDynamic<double> sol, constr;

    engine.SetPMatrix(chronoP);
    engine.SetGMatrix(chronoG);
    engine.SetqMatrix(chrono_q);
    engine.SethMatrix(chrono_h);
    engine.SetNumConstrUnilateralLinear(3);
    engine.SetNumConstrUnilateralConic(std::vector<int_t>(1,4));

    engine.CheckMatrices();

    engine.Run();

    engine.GetSolution(sol);
    engine.GetSlackVariable(constr);

    chrono::GetLog() << sol << "\n";
    chrono::GetLog() << constr << "\n";

    getchar();
    return 0;
}

int main(int argc, char* argv[]) {
    demo_withEngine();
    
}