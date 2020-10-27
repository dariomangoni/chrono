#ifndef CHMODELREDUCTION_H
#define CHMODELREDUCTION_H

// Chrono header files
#include "ChApiModelReduction.h"
#include "chrono/core/ChMatrix.h"

// Mumps header files
#ifdef CHRONO_MUMPS
#include "chrono_mumps/ChMumpsEngine.h"
#endif

// Eigen header files
#include <Eigen/SparseCore>
#include <Eigen/PardisoSupport>

// Spectra header files
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/Util/CompInfo.h>

namespace chrono {
using namespace Spectra;

class ChApiModelReduction ChSymGEigsSolver
{
public:
    ChSymGEigsSolver(const ChSparseMatrix& matA_in,
                     const ChSparseMatrix& matB_in,
                     ChMatrixDynamic<double>& eig_val_out,
                     ChMatrixDynamic<double>& eig_vect_out)
        : matA(matA_in),
        matB(matB_in),
        eig_val(eig_val_out),
        eig_vect(eig_vect_out) {}

    ~ChSymGEigsSolver() {}

    void SetVerbose(bool val) { verbose = val; }

    int computeShift(int requested_eigval, double sigma) const;

    int computeCholesky(int requested_eigval) const;

    int computeRegularInverse(int requested_eigval) const;

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


template <typename Scalar_,
          int UploA = Eigen::Lower,
          int UploB = Eigen::Lower,
          int FlagsA = Eigen::ColMajor,
          int FlagsB = Eigen::ColMajor,
          typename StorageIndexA = int,
          typename StorageIndexB = int>
class SymShiftInvertPardisoLU {
  public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

  private:
    using Index = Eigen::Index;

    // Hypothetical type of the A matrix, either dense or sparse
    using SparseTypeA = Eigen::SparseMatrix<Scalar, FlagsA, StorageIndexA>;

    // Hypothetical type of the B matrix, either dense or sparse
    using SparseTypeB = Eigen::SparseMatrix<Scalar, FlagsB, StorageIndexB>;

    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;

    using ConstGenericMatrixA = const Eigen::Ref<const SparseTypeA>;
    using ConstGenericMatrixB = const Eigen::Ref<const SparseTypeB>;

    ConstGenericMatrixA m_matA;
    ConstGenericMatrixB m_matB;
    const Index m_n;
    Eigen::PardisoLU<SparseTypeA> m_solver;

  public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param A A dense or sparse matrix object, whose type can be `Eigen::Matrix<...>`,
    ///          `Eigen::SparseMatrix<...>`, `Eigen::Map<Eigen::Matrix<...>>`,
    ///          `Eigen::Map<Eigen::SparseMatrix<...>>`, `Eigen::Ref<Eigen::Matrix<...>>`,
    ///          `Eigen::Ref<Eigen::SparseMatrix<...>>`, etc.
    /// \param B A dense or sparse matrix object.
    ///
    SymShiftInvertPardisoLU(ConstGenericMatrixA& A, ConstGenericMatrixB& B) : m_matA(A), m_matB(B), m_n(A.rows()) {

        if (m_n != A.cols() || m_n != B.rows() || m_n != B.cols())
            throw std::invalid_argument("SymShiftInvert: A and B must be square matrices of the same size");
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigma) {
        using SpMat = typename ConstGenericMatrixA::PlainObject;
        SpMat matA = m_matA.template selfadjointView<UploA>();
        SpMat matB = m_matB.template selfadjointView<UploB>();
        SpMat mat = matA - sigma * matB;
        // PardisoLU solver
        m_solver.compute(mat);

        if (!m_solver.info() == Eigen::Success)
            throw std::invalid_argument("SymShiftInvert: factorization failed with the given shift");
    }

    ///
    /// Perform the shift-invert operation \f$y=(A-\sigma B)^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * B) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_solver.solve(x);
    }
};



///
template <typename Scalar_, int Uplo = Eigen::Lower, int Flags = Eigen::ColMajor, typename StorageIndex = int>
class SparseRegularInversePardisoLU {
  public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

  private:
    using Index = Eigen::Index;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar, Flags, StorageIndex>;
    using ConstGenericSparseMatrix = const Eigen::Ref<const SparseMatrix>;

    ConstGenericSparseMatrix m_mat;
    const Index m_n;
    Eigen::PardisoLU<SparseMatrix> m_solver;
    mutable Spectra::CompInfo m_info;

  public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param mat An **Eigen** sparse matrix object, whose type can be
    /// `Eigen::SparseMatrix<Scalar, ...>` or its mapped version
    /// `Eigen::Map<Eigen::SparseMatrix<Scalar, ...> >`.
    ///
    SparseRegularInversePardisoLU(ConstGenericSparseMatrix& mat) : m_mat(mat), m_n(mat.rows()) {
        if (mat.rows() != mat.cols())
            throw std::invalid_argument("SparseRegularInverse: matrix must be square");

        m_solver.compute(mat);
        m_info =
            (m_solver.info() == Eigen::Success) ? Spectra::CompInfo::Successful : Spectra::CompInfo::NumericalIssue;
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Returns the status of the computation.
    /// The full list of enumeration values can be found in \ref Enumerations.
    ///
    Spectra::CompInfo info() const { return m_info; }

    ///
    /// Perform the solving operation \f$y=B^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(B) * x_in
    void solve(const Scalar* x_in, Scalar* y_out) const {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_solver.solve(x);

        m_info = (m_solver.info() == Eigen::Success) ? Spectra::CompInfo::Successful : Spectra::CompInfo::NotConverging;
        if (m_info != Spectra::CompInfo::Successful)
            throw std::runtime_error("SparseRegularInverse: CG solver does not converge");
    }

    ///
    /// Perform the matrix-vector multiplication operation \f$y=Bx\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = B * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
    }
};

#ifdef CHRONO_MUMPS
template <typename Scalar_,
          int UploA = Eigen::Lower,
          int UploB = Eigen::Lower,
          int FlagsA = Eigen::ColMajor,
          int FlagsB = Eigen::ColMajor,
          typename StorageIndexA = int,
          typename StorageIndexB = int>
class SymShiftInvertMumps {
  public:
    ///
    /// Element type of the matrix.
    ///
    using Scalar = Scalar_;

  private:
    using Index = Eigen::Index;

    // Hypothetical type of the A matrix, either dense or sparse
    using SparseTypeA = Eigen::SparseMatrix<Scalar, FlagsA, StorageIndexA>;

    // Hypothetical type of the B matrix, either dense or sparse
    using SparseTypeB = Eigen::SparseMatrix<Scalar, FlagsB, StorageIndexB>;

    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using MapConstVec = Eigen::Map<const Vector>;
    using MapVec = Eigen::Map<Vector>;

    using ConstGenericMatrixA = const Eigen::Ref<const SparseTypeA>;
    using ConstGenericMatrixB = const Eigen::Ref<const SparseTypeB>;

    ConstGenericMatrixA m_matA;
    ConstGenericMatrixB m_matB;
    const Index m_n;
    mutable ChMumpsEngine m_solver;

  public:
    ///
    /// Constructor to create the matrix operation object.
    ///
    /// \param A A dense or sparse matrix object, whose type can be `Eigen::Matrix<...>`,
    ///          `Eigen::SparseMatrix<...>`, `Eigen::Map<Eigen::Matrix<...>>`,
    ///          `Eigen::Map<Eigen::SparseMatrix<...>>`, `Eigen::Ref<Eigen::Matrix<...>>`,
    ///          `Eigen::Ref<Eigen::SparseMatrix<...>>`, etc.
    /// \param B A dense or sparse matrix object.
    ///
    SymShiftInvertMumps(ConstGenericMatrixA& A, ConstGenericMatrixB& B) : m_matA(A), m_matB(B), m_n(A.rows()) {
        m_solver.EnableNullPivotDetection(true);
        m_solver.SetICNTL(14, 150);
        m_solver.SetICNTL(23, 2000);
        m_solver.SetMatrixSymmetry(ChMumpsEngine::mumps_SYM::SYMMETRIC_POSDEF);
        if (m_n != A.cols() || m_n != B.rows() || m_n != B.cols())
            throw std::invalid_argument("SymShiftInvert: A and B must be square matrices of the same size");
    }

    ///
    /// Return the number of rows of the underlying matrix.
    ///
    Index rows() const { return m_n; }
    ///
    /// Return the number of columns of the underlying matrix.
    ///
    Index cols() const { return m_n; }

    ///
    /// Set the real shift \f$\sigma\f$.
    ///
    void set_shift(const Scalar& sigma) {
        using SpMat = typename ConstGenericMatrixA::PlainObject;
        SpMat matA = m_matA.template selfadjointView<UploA>();
        SpMat matB = m_matB.template selfadjointView<UploB>();
        SpMat mat = matA - sigma * matB;
        // Mumps solver
        m_solver.SetMatrix(mat);
        int return_value = m_solver.MumpsCall(ChMumpsEngine::mumps_JOB::ANALYZE_FACTORIZE);

        if (return_value)
            throw std::invalid_argument("SymShiftInvertMumps: factorization failed with the given shift");
    }

    ///
    /// Perform the shift-invert operation \f$y=(A-\sigma B)^{-1}x\f$.
    ///
    /// \param x_in  Pointer to the \f$x\f$ vector.
    /// \param y_out Pointer to the \f$y\f$ vector.
    ///
    // y_out = inv(A - sigma * B) * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const {
        MapConstVec x(x_in, m_n);
        MapVec y(y_out, m_n);
        y = x;
        m_solver.SetRhsVector(y);
        m_solver.MumpsCall(ChMumpsEngine::mumps_JOB::SOLVE);
    }
};
#endif


}

#endif