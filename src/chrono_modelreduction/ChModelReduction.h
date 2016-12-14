#ifndef CHMODELREDUCTION_H
#define CHMODELREDUCTION_H

// Chrono header files
#include "ChApiModelReduction.h"
#include "core/ChSparseMatrix.h"
#include "core/ChMatrixDynamic.h"
// Eigen header files
#include "Eigen/SparseCore"
// Spectra header files
#include <MatOp/SparseSymMatProd.h>
#include "core/ChCSR3Matrix.h"

namespace chrono {
using namespace Spectra;

class ChApiModelReduction ChSymGEigsSolver
{
public:
    ChSymGEigsSolver(const ChCSR3Matrix& matA_in,
                     const ChCSR3Matrix& matB_in,
                     ChMatrixDynamic<double>& eig_val_out,
                     ChMatrixDynamic<double>& eig_vect_out)
        : matA(matA_in),
        matB(matB_in),
        eig_val(eig_val_out),
        eig_vect(eig_vect_out) {}

    ~ChSymGEigsSolver() {}

    int compute(int requested_eigval) const;

protected:

    static Eigen::Map<Eigen::SparseMatrix<double>> getEigenMapSparseMatrix(const ChCSR3Matrix& mat)
    {
        return Eigen::Map<Eigen::SparseMatrix<double>>(mat.GetNumRows(),
                                                       mat.GetNumColumns(),
                                                       mat.GetNNZ(),
                                                       mat.GetCSR_LeadingIndexArray(),
                                                       mat.GetCSR_TrailingIndexArray(),
                                                       mat.GetCSR_ValueArray());
    }


private:
    const ChCSR3Matrix &matA, &matB;
    ChMatrixDynamic<double> &eig_val, &eig_vect;

};


}

#endif