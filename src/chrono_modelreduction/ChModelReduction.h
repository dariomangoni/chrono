#ifndef CHMODELREDUCTION_H
#define CHMODELREDUCTION_H

#include <Eigen/Core>

#include "ChApiModelReduction.h"
#include "core/ChCSR3Matrix.h"
#include <Eigen/src/SparseCore/CompressedStorage.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <MatOp/SparseCholesky.h>

#define EIGEN_SPARSEMATRIX_PLUGIN "chrono_modelreduction/ChModelReductionEigenExtension.h"

class ChApiModelReduction ChCSR3MatrixREDUCTION
{
private:
    //std::shared_ptr<chrono::ChCSR3Matrix> m_mat;
    chrono::ChCSR3Matrix* m_mat;

public:
    ChCSR3MatrixREDUCTION(chrono::ChCSR3Matrix& mat_source) : m_mat(&mat_source) {}
    virtual ~ChCSR3MatrixREDUCTION() {}

    int rows() const { return m_mat->GetNumRows(); }
    int cols() const { return m_mat->GetNumColumns(); }

    // y_out = M * x_in
    void perform_op(double *x_in, double *y_out)
    {
        m_mat = const_cast<chrono::ChCSR3Matrix*>(m_mat->MultiplyVect(x_in, y_out));
    }
};

class ChApiModelReduction ChEigenSparseMatrix :
    public Eigen::SparseMatrix<double, Eigen::RowMajor, int>, public chrono::ChSparseMatrix
{
public:
    ChEigenSparseMatrix()
        : SparseMatrix<double, Eigen::RowMajor, int>(),
        chrono::ChSparseMatrix() {}
    ChEigenSparseMatrix(int rows, int cols)
        : SparseMatrix<double, Eigen::RowMajor, int>(rows, cols),
        chrono::ChSparseMatrix(rows, cols) {}

    // TODO: check consistency of the constructors below
    template<typename OtherDerived>
    ChEigenSparseMatrix(const Eigen::SparseMatrixBase<OtherDerived>& other)
        : Eigen::SparseMatrix<double, Eigen::RowMajor, int>(other),
        chrono::ChSparseMatrix(other.rows(), other.cols()) {}

    //template<typename OtherDerived>
    //ChEigenSparseMatrix& operator=(const Eigen::MatrixBase <OtherDerived>& other)
    //{
    //    this->Eigen::SparseMatrix<double, Eigen::RowMajor, int>::operator=(other);
    //    return *this;
    //}


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
    int* GetCSR_LeadingIndexArray() const override { return outerIndexPtr(); }

    /// Return the column index array in the CSR representation of this matrix.
    int* GetCSR_TrailingIndexArray() const override{ return innerIndexPtr(); }

    /// Return the array of matrix values in the CSR representation of this matrix.
    double* GetCSR_ValueArray() const override { return valuePtr(); }


};

//
////template<typename _Scalar, typename _Index>
////class ChApiModelReduction CompressedStorage_backdoor : public Eigen::CompressedStorage<_Scalar, _Index>
////{
////    CompressedStorage_backdoor() : Eigen::CompressedStorage<_Scalar, _Index>() {};
////    virtual ~CompressedStorage_backdoor()
////    {
////        m_size = 0;
////        m_allocatedSize = 0;
////        m_indices = nullptr;
////        m_values = nullptr;
////    }
////
////    template<typename ChCSR3MatrixIMPORTED>
////    void LoadChCSR3Matrix(ChCSR3MatrixIMPORTED& ch_mat)
////    {
////        ch_mat.Compress();
////
////        m_allocatedSize = ChCSR3MatrixIMPORTED.GetTrailingIndexCapacity();
////        m_size = ChCSR3MatrixIMPORTED.GetTrailingIndexLength();
////        m_indices = ChCSR3MatrixIMPORTED.GetCSR_TrailingIndexArray();
////        m_values = ch_mat.GetCSR_ValueArray();
////
////
////
////    }
////
////};
//
//
int ChApiModelReduction model_reduction_test();
int ChApiModelReduction model_reduction_test2();
int ChApiModelReduction model_reduction_usermultiplication();
int ChApiModelReduction model_reduction_generalized();
int ChApiModelReduction model_reduction_generalized2();

#endif