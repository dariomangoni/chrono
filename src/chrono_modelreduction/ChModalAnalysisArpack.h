#ifndef CHARPACKSOLVER_H
#define CHARPACKSOLVER_H

// Chrono header files
#include "ChApiModelReduction.h"
#include "chrono/core/ChMatrix.h"
#include "arpack.hpp"

// Eigen header files
#include <Eigen/SparseCore>
#include <Eigen/PardisoSupport>


namespace chrono {

class ChApiModelReduction ChArpackSolver
{
public:
    ChArpackSolver(const ChSparseMatrix& matA_in,
                     const ChSparseMatrix& matB_in,
                     ChMatrixDynamic<double>& eig_val_out,
                     ChMatrixDynamic<double>& eig_vect_out)
        : matA(matA_in),
        matB(matB_in),
        eig_val(eig_val_out),
        eig_vect(eig_vect_out) {}

    ~ChArpackSolver() {}

    void SetVerbose(bool val) { verbose = val; }

    int compute(int requested_eigval, double sigmar, double sigmai) const;

    /// TODO: this functions returns a sparse matrix in column major order as needed by Eigen;
    /// in fact SimplicialLLDL used internally by Spectra fails with row major matrices (as ChSparseMatrix)
    static Eigen::Map<Eigen::SparseMatrix<double>> getEigenMapSparseMatrix(const ChSparseMatrix& mat)
    {

        return Eigen::Map<Eigen::SparseMatrix<double>>(const_cast<ChSparseMatrix&>(mat).rows(),
                                                       const_cast<ChSparseMatrix&>(mat).cols(),
                                                       const_cast<ChSparseMatrix&>(mat).nonZeros(),
                                                       const_cast<ChSparseMatrix&>(mat).outerIndexPtr(),
                                                       const_cast<ChSparseMatrix&>(mat).innerIndexPtr(),
                                                       const_cast<ChSparseMatrix&>(mat).valuePtr());
    }

private:

    const ChSparseMatrix &matA, &matB;
    ChMatrixDynamic<double> &eig_val, &eig_vect;
    bool verbose = false;

};


}

#endif