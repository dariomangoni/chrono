#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <iostream>
using namespace Spectra;
int main()
{
    // We are going to solve the generalized eigenvalue problem A * x = lambda * B * x
    const int n = 100;
    // Define the A matrix
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd A = M + M.transpose();
    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    Eigen::SparseMatrix<double> B(n, n);
    B.reserve(Eigen::VectorXi::Constant(n, 3));
    for(int i = 0; i < n; i++)
    {
        B.insert(i, i) = 2.0;
        if(i > 0)
            B.insert(i - 1, i) = 1.0;
        if(i < n - 1)
            B.insert(i + 1, i) = 1.0;
    }
    // Construct matrix operation object using the wrapper classes
    DenseSymMatProd<double> op(A);
    SparseCholesky<double>  Bop(B);
    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    SymGEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY>
        geigs(&op, &Bop, 3, 6);
    // Initialize and compute
    geigs.init();
    int nconv = geigs.compute();
    // Retrieve results
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evecs;
    if(geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }
    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;
    // Verify results using the generalized eigen solver in Eigen
    Eigen::MatrixXd Bdense = B;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Bdense);
    std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    return 0;
}