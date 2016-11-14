#include <Eigen/Core>
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseSymMatProd.h>
#include <iostream>
#include "ChModelReduction.h"
#include <core/ChCSR3Matrix.h>
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


// M = diag(1, 2, ..., 10)
class ChCSR3MatrixREDUCTION
{
private:
    //std::shared_ptr<chrono::ChCSR3Matrix> m_mat;
    chrono::ChCSR3Matrix* m_mat;

public:
    ChCSR3MatrixREDUCTION(chrono::ChCSR3Matrix& mat_source) : m_mat(&mat_source) {}
    ~ChCSR3MatrixREDUCTION(){}

    int rows() const { return m_mat->GetNumRows(); }
    int cols() const { return m_mat->GetNumColumns(); }

    // y_out = M * x_in
    void perform_op(double *x_in, double *y_out)
    {
        m_mat = const_cast<chrono::ChCSR3Matrix*>(m_mat->MultiplyVect(x_in, y_out));
    }
};

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