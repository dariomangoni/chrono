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

#include <bitset>

#include "chrono_cvxopt/ChCvxoptConeQpEngine.h"

namespace chrono {

ChCvxoptConeQpEngine::ChCvxoptConeQpEngine() {
    Py_Initialize();

    if (import_cvxopt() < 0) {
        throw ChException("Chrono::cvxopt: cannot import cvxopt library.");
    }

    /* import cvxopt.solvers */
    solverModule = PyImport_ImportModule("cvxopt.solvers");
    if (!solverModule) {
        throw ChException("Chrono::cvxopt: cannot import cvxopt.solver module.");
    }

    /* get reference to solvers.solvelp */
    programSolver = PyObject_GetAttrString(solverModule, "coneqp");
    if (!programSolver) {
        Py_DECREF(solverModule);
        throw ChException("Chrono::cvxopt: cannot reference cvxopt.solvers.coneqp");
    }

    // objects are created, but their size are zero
    // they are just wrappers of the real matrices
    P = SpMatrix_New(0, 0, 0, DOUBLE);
    G = SpMatrix_New(0, 0, 0, DOUBLE);
    A = SpMatrix_New(0, 0, 0, DOUBLE);
    q = Matrix_New(0, 1, DOUBLE);
    h = Matrix_New(0, 1, DOUBLE);
    b = Matrix_New(0, 1, DOUBLE);

    // pointers are increased so that no dealloc will be called on the underlying structures
    Py_INCREF(P);
    Py_INCREF(G);
    Py_INCREF(A);
    Py_INCREF(q);
    Py_INCREF(h);
    Py_INCREF(b);

    PyDict_SetItemString(dims, "s", PyList_New(0));  // semidefinite cones; not implemented
    
}

void ChCvxoptConeQpEngine::SetPMatrix(ChCSMatrixIndexType<int_t>& Pmat) {
    auto Pccs_ptr = reinterpret_cast<spmatrix*>(P)->obj;
    Pccs_ptr->colptr = Pmat.GetCS_LeadingIndexArray();
    Pccs_ptr->rowind = Pmat.GetCS_TrailingIndexArray();
    Pccs_ptr->values = Pmat.GetCS_ValueArray();
    Pccs_ptr->nrows = Pmat.GetNumRows();
    Pccs_ptr->ncols = Pmat.GetNumColumns();
}

void ChCvxoptConeQpEngine::SetGMatrix(ChCSMatrixIndexType<int_t>& Gmat) {
    auto Gccs_ptr = reinterpret_cast<spmatrix*>(G)->obj;
    Gccs_ptr->colptr = Gmat.GetCS_LeadingIndexArray();
    Gccs_ptr->rowind = Gmat.GetCS_TrailingIndexArray();
    Gccs_ptr->values = Gmat.GetCS_ValueArray();
    Gccs_ptr->nrows = Gmat.GetNumRows();
    Gccs_ptr->ncols = Gmat.GetNumColumns();
}

void ChCvxoptConeQpEngine::SetAMatrix(ChCSMatrixIndexType<int_t>& Amat) {
    auto Accs_ptr = reinterpret_cast<spmatrix*>(A)->obj;
    Accs_ptr->colptr = Amat.GetCS_LeadingIndexArray();
    Accs_ptr->rowind = Amat.GetCS_TrailingIndexArray();
    Accs_ptr->values = Amat.GetCS_ValueArray();
    Accs_ptr->nrows = Amat.GetNumRows();
    Accs_ptr->ncols = Amat.GetNumColumns();
}

void ChCvxoptConeQpEngine::SetqMatrix(ChMatrix<>& qmat) {
    auto q_ccs_ptr = reinterpret_cast<matrix*>(q);
    q_ccs_ptr->nrows = qmat.GetRows();
    q_ccs_ptr->ncols = qmat.GetColumns();
    q_ccs_ptr->buffer = qmat.GetAddress();
}

void ChCvxoptConeQpEngine::SethMatrix(ChMatrix<>& hmat) {
    auto h_ccs_ptr = reinterpret_cast<matrix*>(h);
    h_ccs_ptr->nrows = hmat.GetRows();
    h_ccs_ptr->ncols = hmat.GetColumns();
    h_ccs_ptr->buffer = hmat.GetAddress();
}

void ChCvxoptConeQpEngine::SetbMatrix(ChMatrix<>& bmat) {
    auto b_ccs_ptr = reinterpret_cast<matrix*>(b);
    b_ccs_ptr->nrows = bmat.GetRows();
    b_ccs_ptr->ncols = bmat.GetColumns();
    b_ccs_ptr->buffer = bmat.GetAddress();
}

void ChCvxoptConeQpEngine::SetNumConstrUnilateralLinear(int_t constr_unilin_dim) {
    PyDict_SetItemString(dims, "l", PyInt_FromSsize_t(constr_unilin_dim));
}

void ChCvxoptConeQpEngine::SetNumConstrUnilateralConic(const std::vector<Py_ssize_t>& unicon_dims) {
    PyObject* constr_unicone_dims_list = PyList_New(unicon_dims.size());
    for (auto cone_set = 0; cone_set < unicon_dims.size(); ++cone_set) {
        PyList_SetItem(constr_unicone_dims_list, cone_set, PyInt_FromSsize_t(unicon_dims.at(cone_set)));
    }

    PyDict_SetItemString(dims, "q", constr_unicone_dims_list);
}

void ChCvxoptConeQpEngine::GetSolution(ChMatrix<>& x) const {
    if (!getDenseMatFromSolution(x, "x"))
        throw ChException("chrono::cvxopt: the solution vector does not exist");
}

void ChCvxoptConeQpEngine::GetSlackVariable(ChMatrix<>& s) const {
    getDenseMatFromSolution(s, "s");
}

void ChCvxoptConeQpEngine::GetLagrangianMultBilateral(ChMatrix<>& y) const {
    getDenseMatFromSolution(y, "y");
}

void ChCvxoptConeQpEngine::GetLagrangianMultUnilateral(ChMatrix<>& z) const {
    getDenseMatFromSolution(z, "z");
}

bool ChCvxoptConeQpEngine::IsSolutionOptimal() const {
    return sol && strcmp(PyString_AsString(PyDict_GetItemString(sol, "status")), "optimal") == 0 ? true : false;
}

void ChCvxoptConeQpEngine::CheckMatrices() const {
    std::cout << "P: " << (SpMatrix_Check(P) ? "OK" : "ERROR") << std::endl;
    std::cout << "G: " << (SpMatrix_Check(G) ? "OK" : "ERROR") << std::endl;
    std::cout << "A: " << (SpMatrix_Check(A) ? "OK" : "ERROR") << std::endl;
    std::cout << "q: " << (Matrix_Check(q) ? "OK" : "ERROR") << std::endl;
    std::cout << "h: " << (Matrix_Check(h) ? "OK" : "ERROR") << std::endl;
    std::cout << "b: " << (Matrix_Check(b) ? "OK" : "ERROR") << std::endl;
}

void ChCvxoptConeQpEngine::Run() {
    // inform CVXOPT about the dimensions of the constraints

    /* pack matrices into an argument tuple - references are stolen*/
    /* since references have been stolen they have to be increased to avoid garbage collection*/
    // TODO: check if this is true
    Py_INCREF(P);
    Py_INCREF(G);
    Py_INCREF(q);
    Py_INCREF(h);
    Py_INCREF(dims);

    PyObject* pArgs = PyTuple_New(SP_NROWS(A) > 0 && MAT_NROWS(b) > 0 ? 7 : 5);
    PyTuple_SetItem(pArgs, 0, reinterpret_cast<PyObject*>(P));
    PyTuple_SetItem(pArgs, 1, reinterpret_cast<PyObject*>(q));
    PyTuple_SetItem(pArgs, 2, reinterpret_cast<PyObject*>(G));
    PyTuple_SetItem(pArgs, 3, reinterpret_cast<PyObject*>(h));
    PyTuple_SetItem(pArgs, 4, dims);

    if (SP_NROWS(A) > 0 && MAT_NROWS(b) > 0) {
        Py_INCREF(A);
        Py_INCREF(b);
        PyTuple_SetItem(pArgs, 5, reinterpret_cast<PyObject*>(A));
        PyTuple_SetItem(pArgs, 6, reinterpret_cast<PyObject*>(b));
    }

    sol = PyObject_CallObject(programSolver, pArgs);
    if (!sol) {
        GetPythonError(true);
        Py_DECREF(solverModule);
        Py_DECREF(programSolver);
    }

    Py_DECREF(pArgs);
}

void ChCvxoptConeQpEngine::GetPythonError(bool print_fulltrace) {
    auto err = PyErr_Occurred();
    if (err != nullptr) {
        std::string full_backtrace;
        PyObject *ptype, *pvalue, *ptraceback;
        PyObject *pystr, *module_name, *pyth_module, *pyth_func;

        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
        pystr = PyObject_Str(pvalue);
        std::string error_description = PyString_AsString(pystr);

        /* See if we can get a full traceback */
        module_name = PyString_FromString("traceback");
        pyth_module = PyImport_Import(module_name);
        Py_DECREF(module_name);

        if (pyth_module != nullptr) {
            pyth_func = PyObject_GetAttrString(pyth_module, "format_exception");
            if (pyth_func && PyCallable_Check(pyth_func)) {
                PyObject* pyth_val;

                pyth_val = PyObject_CallFunctionObjArgs(pyth_func, ptype, pvalue, ptraceback, NULL);

                pystr = PyObject_Str(pyth_val);
                full_backtrace = PyString_AsString(pystr);
                Py_DECREF(pyth_val);
            }
        }

        std::cout << "Error: " << error_description << std::endl;
        if (print_fulltrace && full_backtrace.length() > 0)
            std::cout << "Full Traceback: \n" << full_backtrace << std::endl;

        Py_DECREF(ptype);
        Py_DECREF(pvalue);
        Py_DECREF(ptraceback);
        Py_DECREF(pystr);
        Py_DECREF(pyth_func);
        Py_DECREF(pyth_module);

        throw ChException("chrono::cvxopt: error during cvxopt coneqp run.");
    }
}

bool ChCvxoptConeQpEngine::getDenseMatFromSolution(ChMatrix<>& matout, std::string key) const {
    PyObject* mat_cvxopt = PyDict_GetItemString(sol, key.c_str());

    if (mat_cvxopt) {
        matout.Reset(MAT_NROWS(mat_cvxopt), MAT_NCOLS(mat_cvxopt));

        for (auto row_sel = 0; row_sel < MAT_NROWS(mat_cvxopt); ++row_sel)
            for (auto col_sel = 0; col_sel < MAT_NCOLS(mat_cvxopt); ++col_sel)
                matout.SetElement(row_sel, col_sel, MAT_BUFD(mat_cvxopt)[row_sel + col_sel * MAT_NROWS(mat_cvxopt)]);
        return true;
    } else {
        matout.Reset(0, 0);
        return false;
    }
}

ChCvxoptConeQpEngine::~ChCvxoptConeQpEngine() {
    Py_DECREF(solverModule);
    Py_DECREF(programSolver);
    Py_XDECREF(sol);

    // reference counters to P, G, etc... are still +1
    // however, before decreasing them they should be resized to 0,0 dimension

    Py_Finalize();
}
}  // namespace chrono
