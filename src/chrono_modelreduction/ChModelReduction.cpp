
#include "ChModelReduction.h"
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>

#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

namespace chrono
{
    int chrono::ChSymGEigsSolver::compute(int requested_eigval) const
    {
        int dim = matA.rows();
        assert(dim > 0 && "Matrix dimensions are not valid");
        int convergence_speed = static_cast<int>(std::max(std::min((dim + 4*requested_eigval) / 2, dim), requested_eigval) );

        // Construct matrix operation object using the wrapper classes
        assert(matA.isCompressed() && "ChSymGEigsSolver: matA is not in a standard CSR format");
        assert(matB.isCompressed() && "ChSymGEigsSolver: matB is not in a standard CSR format");
        
        int m = convergence_speed;
        double sigma = 100;

        using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
        using BOpType = SparseSymMatProd<double>;

        OpType op(getEigenMapSparseMatrix(matA), getEigenMapSparseMatrix(matB));
        BOpType Bop(getEigenMapSparseMatrix(matB));
        SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> geigs(op, Bop, requested_eigval, m, sigma);


        //SparseSymMatProd<double> op(getEigenMapSparseMatrix(matA));
        //SparseCholesky<double> Bop(getEigenMapSparseMatrix(matB));

        ////SparseSymMatProd<double> op(matA);
        ////SparseCholesky<double> Bop(matB);


        //// Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
        //const SELECT_EIGENVALUE SelectionRule = SMALLEST_MAGN;
        //SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, requested_eigval, convergence_speed);

        // Initialize and compute
        geigs.init();
        int nconv = geigs.compute(SortRule::LargestMagn);

        // Retrieve results
        eig_val.resize(requested_eigval, 1);
        eig_vect.resize(dim, requested_eigval);
        Eigen::Map<Eigen::VectorXd> eig_val_map(eig_val.data(), requested_eigval, 1); //TODO: use default Chrono Vectors
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eig_vect_map(eig_vect.data(), dim, requested_eigval);

        if (geigs.info() == CompInfo::Successful)
        {
            eig_val_map = geigs.eigenvalues();
            eig_vect_map = geigs.eigenvectors();
        }
        
        return geigs.info() == CompInfo::Successful ? 0:1;
    }
}


