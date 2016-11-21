#ifndef CHMODELREDUCTION_H
#define CHMODELREDUCTION_H

#include "ChApiModelReduction.h"
#include "core/ChCSR3Matrix.h"
#include <Eigen/src/SparseCore/CompressedStorage.h>

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

//template<typename _Scalar, typename _Index>
//class ChApiModelReduction CompressedStorage_backdoor : public Eigen::CompressedStorage<_Scalar, _Index>
//{
//    CompressedStorage_backdoor() : Eigen::CompressedStorage<_Scalar, _Index>() {};
//    virtual ~CompressedStorage_backdoor()
//    {
//        m_size = 0;
//        m_allocatedSize = 0;
//        m_indices = nullptr;
//        m_values = nullptr;
//    }
//
//    template<typename ChCSR3MatrixIMPORTED>
//    void LoadChCSR3Matrix(ChCSR3MatrixIMPORTED& ch_mat)
//    {
//        ch_mat.Compress();
//
//        m_allocatedSize = ChCSR3MatrixIMPORTED.GetTrailingIndexCapacity();
//        m_size = ChCSR3MatrixIMPORTED.GetTrailingIndexLength();
//        m_indices = ChCSR3MatrixIMPORTED.GetCSR_TrailingIndexArray();
//        m_values = ch_mat.GetCSR_ValueArray();
//
//
//
//    }
//
//};


int ChApiModelReduction model_reduction_test();
int ChApiModelReduction model_reduction_test2();
int ChApiModelReduction model_reduction_usermultiplication();
int ChApiModelReduction model_reduction_generalized();
int ChApiModelReduction model_reduction_generalized2();

#endif