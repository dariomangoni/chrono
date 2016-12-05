
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseSymMatProd.h>
#include <iostream>
#include "ChModelReduction.h"
#include <GenEigsSolver.h>


#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random> // Requires C++ 11

#include <SymGEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/DenseCholesky.h>
#include <MatOp/SparseSymMatProd.h>
#include <MatOp/SparseCholesky.h>
#include "core/ChCSR3Matrix.h"
using namespace Spectra;



int model_reduction_test()
{
    // We are going to calculate the eigenvalues of M
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
    Eigen::MatrixXd M = A + A.transpose();
    // Construct matrix operation object using the wrapper class DenseGenMatProd
    DenseSymMatProd<double> op(M);
    // Construct eigen solver object, requesting the largest three eigenvalues
    SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 3, 6);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(1000, 1e-12, LARGEST_ALGE);
    // Retrieve results
    Eigen::VectorXd evalues;
    if (eigs.info() == SUCCESSFUL)
        evalues = eigs.eigenvalues();
    std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    return 0;
}

int model_reduction_test2()
{
    // We are going to calculate the eigenvalues of M
    Eigen::MatrixXd M(5,5);
    M.resize(5, 5);
    for (auto row_sel = 0; row_sel<M.rows(); ++row_sel)
    {
        for (auto col_sel = 0; col_sel<M.cols(); ++col_sel)
        {
            M.coeffRef(row_sel, col_sel) = 0;
        }
    }
    M.coeffRef(0, 0) = 1;
    M.coeffRef(1, 1) = 2;
    M.coeffRef(2, 2) = 3;
    M.coeffRef(3, 3) = 4;
    M.coeffRef(4, 4) = 5;
    // Construct matrix operation object using the wrapper class DenseGenMatProd
    DenseSymMatProd<double> op(M);
    // Construct eigen solver object, requesting the largest three eigenvalues
    SymEigsSolver< double, SMALLEST_MAGN, DenseSymMatProd<double> > eigs(&op, 3, 4);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(1000, 1e-14, SMALLEST_MAGN);
    // Retrieve results
    Eigen::VectorXd evalues;
    if (eigs.info() == SUCCESSFUL)
        evalues = eigs.eigenvalues();
    std::cout << "Eigenvalues found:\n" << evalues << std::endl;

    std::cout << "M is:" << std::endl;
    for (auto row_sel = 0; row_sel<M.rows(); ++row_sel)
    {
        for (auto col_sel = 0; col_sel<M.cols(); ++col_sel)
        {
            std::cout << M.coeff(row_sel, col_sel) << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}




int model_reduction_usermultiplication()
{
    const int n = 5;
    chrono::ChCSR3Matrix mymat(n, n);

    mymat.SetElement(0, 0, 1);
    mymat.SetElement(1, 1, 2);
    mymat.SetElement(2, 2, 3);
    mymat.SetElement(3, 3, 4);
    mymat.SetElement(4, 4, 5);
    mymat.Compress();

    ChCSR3MatrixREDUCTION mymat_red(mymat);

    std::cout << "mymat is:" << std::endl;
    for (auto row_sel = 0; row_sel<mymat.GetNumRows(); ++row_sel)
    {
        for (auto col_sel = 0; col_sel<mymat.GetNumColumns(); ++col_sel)
        {
            std::cout << mymat.GetElement(row_sel, col_sel) << " ";
        }
        std::cout << std::endl;
    }

    SymEigsSolver<double, SMALLEST_MAGN, ChCSR3MatrixREDUCTION> eigs(&mymat_red, 1, 3);
    eigs.init();
    eigs.compute();
    if (eigs.info() == SUCCESSFUL)
    {
        Eigen::VectorXd evalues = eigs.eigenvalues();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    }
    else
        std::cout << "Something went wrong\n" << std::endl;

    return 0;
}



int model_reduction_onlyeigen()
{
    const int n = 5;
    chrono::ChCSR3Matrix mymat(n, n);

    mymat.SetElement(0, 0, 1);
    mymat.SetElement(1, 1, 2);
    mymat.SetElement(2, 2, 3);
    mymat.SetElement(3, 3, 4);
    mymat.SetElement(4, 4, 5);
    mymat.Compress();

    ChCSR3MatrixREDUCTION mymat_red(mymat);

    std::cout << "mymat is:" << std::endl;
    for (auto row_sel = 0; row_sel<mymat.GetNumRows(); ++row_sel)
    {
        for (auto col_sel = 0; col_sel<mymat.GetNumColumns(); ++col_sel)
        {
            std::cout << mymat.GetElement(row_sel, col_sel) << " ";
        }
        std::cout << std::endl;
    }

    SymEigsSolver<double, SMALLEST_MAGN, ChCSR3MatrixREDUCTION> eigs(&mymat_red, 1, 3);
    eigs.init();
    eigs.compute();
    if (eigs.info() == SUCCESSFUL)
    {
        Eigen::VectorXd evalues = eigs.eigenvalues();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    }
    else
        std::cout << "Something went wrong\n" << std::endl;

    return 0;
}

// // -------------------------------------------------------------------------


using namespace Spectra;
    
#define CATCH_CONFIG_MAIN
#include "../test/catch.hpp"

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::SparseMatrix<double> SpMatrix;

// Traits to obtain operation type from matrix type
template <typename MatType>
struct OpTypeTrait
{
    typedef DenseSymMatProd<double> OpType;
};

template <>
struct OpTypeTrait<SpMatrix>
{
    typedef SparseSymMatProd<double> OpType;
};

template <typename MatType>
struct BOpTypeTrait
{
    typedef DenseCholesky<double> OpType;
};

template <>
struct BOpTypeTrait<SpMatrix>
{
    typedef SparseCholesky<double> OpType;
};



// Generate data for testing
void gen_dense_data(int n, Matrix& A, Matrix& B)
{
    Matrix M = Eigen::MatrixXd::Random(n, n);
    A = M + M.transpose();
    B = M.transpose() * M;
    // To make sure B is positive definite
    B.diagonal() += Eigen::VectorXd::Random(n).cwiseAbs();
}


int model_reduction_generalized()
{

    //TEST_CASE("Generalized eigensolver of symmetric real matrix [10x10]", "[geigs_sym]")
    
    std::srand(123);
    bool allow_fail = true;
    const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;

    Matrix A, B;
    gen_dense_data(10, A, B);
    int k = 3;
    int m = 6;

    typedef typename OpTypeTrait<Matrix>::OpType OpType;
    typedef typename BOpTypeTrait<Matrix>::OpType BOpType;
    OpType op(A);
    BOpType Bop(B);
    SymGEigsSolver<double, SelectionRule, OpType, BOpType, GEIGS_CHOLESKY> eigs(&op, &Bop, k, m);
    eigs.init();
    int nconv = eigs.compute(100); // maxit = 100 to reduce running time for failed cases
    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    if (allow_fail)
    {
        if (eigs.info() != SUCCESSFUL)
        {
            std::cout << "FAILED on this test";
            std::cout << "nconv = " << nconv << std::endl;
            std::cout << "niter = " << niter << std::endl;
            std::cout << "nops  = " << nops << std::endl;
            return 1;
        }
    }
    else {
        std::cout << "nconv = " << nconv;
        std::cout << "niter = " << niter;
        std::cout << "nops  = " << nops;
        assert(eigs.info() == SUCCESSFUL);
    }

    Vector evals = eigs.eigenvalues();
    Matrix evecs = eigs.eigenvectors();

    Matrix resid = A.template selfadjointView<Eigen::Lower>() * evecs -
        B.template selfadjointView<Eigen::Lower>() * evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    std::cout << "||AU - BUD||_inf = " << err;
    assert(err == Approx(0.0));


    return 0;
}


// ------------------------------------------
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <SymGEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseCholesky.h>
#include <iostream>
using namespace Spectra;
int model_reduction_generalized2()
{
    // We are going to solve the generalized eigenvalue problem A * x = lambda * B * x
    const int n = 100;
    // Define the A matrix
    Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd A = M + M.transpose();
    // Define the B matrix, a band matrix with 2 on the diagonal and 1 on the subdiagonals
    Eigen::SparseMatrix<double> B(n, n);
    B.reserve(Eigen::VectorXi::Constant(n, 3));
    for (int i = 0; i < n; i++)
    {
        B.insert(i, i) = 2.0;
        if (i > 0)
            B.insert(i - 1, i) = 1.0;
        if (i < n - 1)
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
    if (geigs.info() == SUCCESSFUL)
    {
        evalues = geigs.eigenvalues();
        evecs = geigs.eigenvectors();
    }
    std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
    std::cout << "Generalized eigenvectors found:\n" << evecs.topRows(10) << std::endl;

    //// Verify results using the generalized eigen solver in Eigen
    //Eigen::MatrixXd Bdense = B;
    //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, Bdense);
    //std::cout << "Generalized eigenvalues:\n" << es.eigenvalues().tail(3) << std::endl;
    //std::cout << "Generalized eigenvectors:\n" << es.eigenvectors().rightCols(3).topRows(10) << std::endl;
    //return 0;
}

