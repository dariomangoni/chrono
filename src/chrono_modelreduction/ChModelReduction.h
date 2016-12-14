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

    int compute(int requested_eigval);

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


    /* *************************************************************************
     * Old stuff bin
     *
     */

class ChEigenSparseMatrix :
    public Eigen::SparseMatrix<double>, public chrono::ChSparseMatrix
{
public:
    ChEigenSparseMatrix()
        : SparseMatrix<double>(),
        chrono::ChSparseMatrix() {}
    ChEigenSparseMatrix(int rows, int cols)
        : SparseMatrix<double>(rows, cols),
        chrono::ChSparseMatrix(rows, cols) {}

    // TODO: check consistency of the constructors below
    template<typename OtherDerived>
    ChEigenSparseMatrix(const Eigen::SparseMatrixBase<OtherDerived>& other)
        : Eigen::SparseMatrix<double>(other),
        chrono::ChSparseMatrix(other.rows(), other.cols()) {}

    template<typename OtherDerived>
    ChEigenSparseMatrix& operator=(const Eigen::MatrixBase <OtherDerived>& other)
    {
        this->Eigen::SparseMatrix<double>::operator=(other);
        return *this;
    }


    // TODO:~SparseMatrix is not virtual!! Memory leaks
    virtual ~ChEigenSparseMatrix() { }

    // ChSparseMatrix operations
    void SetElement(int insrow, int inscol, double insval, bool overwrite = true) override {
        overwrite ? coeffRef(insrow, inscol) = insval : coeffRef(insrow, inscol) += insval;
    }

    double GetElement(int row, int col) const override {
        return coeff(row, col);
    }

    void Reset(int row, int col, int nonzeros = 0) override
    {
        resize(row, col);
        reserve(nonzeros); //TODO: may not be correct
    }

    bool Resize(int nrows, int ncols, int nonzeros = 0) override {
        Reset(nrows, ncols, nonzeros);
        return true;
    }

    bool Compress() override
    {
        makeCompressed();
        return true;
    }

    // TODO: check if outer and inner are not the other way around :-)
    /// Return the row index array in the CSR representation of this matrix.
    int* GetCSR_LeadingIndexArray() const override { return const_cast<int*>(outerIndexPtr()); }

    /// Return the column index array in the CSR representation of this matrix.
    int* GetCSR_TrailingIndexArray() const override { return const_cast<int*>(innerIndexPtr()); }

    /// Return the array of matrix values in the CSR representation of this matrix.
    double* GetCSR_ValueArray() const override { return const_cast<double*>(valuePtr()); }

};


class ChEigenSparseMatrixWrapper : public chrono::ChSparseMatrix
{
private:
    /// matrix in Eigen::SparseMatrix format; it is ColMajor, but since the matrix is symmetric, the CS arrays are the same.
    Eigen::SparseMatrix<double> eig_mat;
public:
    ChEigenSparseMatrixWrapper() :
        chrono::ChSparseMatrix(), eig_mat() {}
    ChEigenSparseMatrixWrapper(int rows, int cols) :
        chrono::ChSparseMatrix(rows, cols), eig_mat(rows, cols) {}

    virtual ~ChEigenSparseMatrixWrapper() { }

    // ChSparseMatrix operations
    void SetElement(int insrow, int inscol, double insval, bool overwrite = true) override {
        overwrite ? eig_mat.coeffRef(insrow, inscol) = insval : eig_mat.coeffRef(insrow, inscol) += insval;
    }

    double GetElement(int row, int col) const override {
        return eig_mat.coeff(row, col);
    }

    void Reset(int row, int col, int nonzeros = 0) override
    {
        eig_mat.resize(row, col);
        eig_mat.reserve(nonzeros); //TODO: may not be correct
    }

    bool Resize(int nrows, int ncols, int nonzeros = 0) override {
        Reset(nrows, ncols, nonzeros);
        return true;
    }

    bool Compress() override
    {
        eig_mat.makeCompressed();
        return true;
    }

    // Eigen operations
    Eigen::SparseMatrix<double>& GetInternalEigenMatrix() { return eig_mat; }

    // TODO: check if outer and inner are not the other way around :-)
    /// Return the row index array in the CSR representation of this matrix.
    int* GetCSR_LeadingIndexArray() const override { return const_cast<int*>(eig_mat.outerIndexPtr()); }

    /// Return the column index array in the CSR representation of this matrix.
    int* GetCSR_TrailingIndexArray() const override { return const_cast<int*>(eig_mat.innerIndexPtr()); }

    /// Return the array of matrix values in the CSR representation of this matrix.
    double* GetCSR_ValueArray() const override { return const_cast<double*>(eig_mat.valuePtr()); }


};



//class ChApiModelReduction ChSymGEigsSolver
//{
//public:
//    ChSymGEigsSolver(ChEigenSparseMatrixWrapper& matA_in,
//                     ChEigenSparseMatrixWrapper& matB_in,
//                     ChMatrixDynamic<double>& eig_val_out,
//                     ChMatrixDynamic<double>& eig_vect_out)
//        : matK(matA_in),
//        matM(matB_in),
//        eig_val(eig_val_out),
//        eig_vect(eig_vect_out) {}
//
//    ~ChSymGEigsSolver() {}
//
//    void compute(int requested_eigval)
//    {
//        int convergence_speed = std::ceil((matK.GetNumRows() + requested_eigval) / 2);
//        // Construct matrix operation object using the wrapper classes
//        SparseSymMatProd<double> op(matK.GetInternalEigenMatrix());
//        SparseCholesky<double>  Bop(matM.GetInternalEigenMatrix());
//
//        // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
//        const SELECT_EIGENVALUE SelectionRule = LARGEST_ALGE;
//        SymGEigsSolver<double, SelectionRule, SparseSymMatProd<double>, SparseCholesky<double>, GEIGS_CHOLESKY> geigs(&op, &Bop, requested_eigval, convergence_speed);
//
//        // Initialize and compute
//        geigs.init();
//        int nconv = geigs.compute();
//
//        // Retrieve results
//        eig_val.Resize(requested_eigval, 1);
//        eig_vect.Resize(matK.GetNumRows(), requested_eigval);
//        Eigen::Map<Eigen::VectorXd> eig_val_map(eig_val.GetAddress(), requested_eigval, 1);
//        Eigen::Map<Eigen::MatrixXd> eig_vect_map(eig_val.GetAddress(), matK.GetNumRows(), requested_eigval);
//
//        if (geigs.info() == SUCCESSFUL)
//        {
//            eig_val_map = geigs.eigenvalues();
//            eig_vect_map = geigs.eigenvectors();
//        }
//
//        std::cout << "Generalized eigenvalues found:\n" << eig_val_map << std::endl;
//        std::cout << "Generalized eigenvectors found:\n" << eig_vect_map.topRows(matK.GetNumRows()) << std::endl;
//    }
//
//private:
//    ChEigenSparseMatrixWrapper& matK, matM;
//    ChMatrixDynamic<double>& eig_val, eig_vect;
//
//};





//int ChApiModelReduction model_reduction_test();
//int ChApiModelReduction model_reduction_test2();
//int ChApiModelReduction model_reduction_usermultiplication();
//int ChApiModelReduction model_reduction_generalized();
//int ChApiModelReduction model_reduction_generalized2();

}

#endif