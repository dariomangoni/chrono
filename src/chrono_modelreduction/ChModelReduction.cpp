
#include "ChModelReduction.h"
#include <SymGEigsSolver.h>
#include <MatOp/SparseCholesky.h>

namespace chrono
{
    int chrono::ChSymGEigsSolver::compute(int requested_eigval)
    {
        int dim = matA.GetNumRows();
        assert(dim > 0 && "Matrix dimensions are not valid");
        int convergence_speed = static_cast<int>(std::max(std::min((dim + 4*requested_eigval) / 2, dim), requested_eigval) );

        // Construct matrix operation object using the wrapper classes
        assert(matA.IsCompressed() && "ChSymGEigsSolver: matA is not in a standard CSR format");
        assert(matB.IsCompressed() && "ChSymGEigsSolver: matB is not in a standard CSR format");
        SparseSymMatProd<double> op(getEigenMapSparseMatrix(matA));
        SparseCholesky<double> Bop(getEigenMapSparseMatrix(matB));

        matA.VerifyMatrix();
        matB.VerifyMatrix();

        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        const SELECT_EIGENVALUE SelectionRule = SMALLEST_MAGN;
        SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, requested_eigval, convergence_speed);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute();

        // Retrieve results
        eig_val.Reset(requested_eigval, 1);
        eig_vect.Reset(dim, requested_eigval);
        Eigen::Map<Eigen::VectorXd> eig_val_map(eig_val.GetAddress(), requested_eigval, 1);
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eig_vect_map(eig_vect.GetAddress(), dim, requested_eigval);

        if (geigs.info() == SUCCESSFUL)
        {
            eig_val_map = geigs.eigenvalues();
            eig_vect_map = geigs.eigenvectors();
        }
        
        return geigs.info();
    }
}


