// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================

#include <numeric>
#include <iomanip>

#include "chrono_modal/ChModalSolverUndamped.h"
#include "chrono_modal/ChKrylovSchurEig.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/solver/ChDirectSolverLScomplex.h"
#include "chrono/utils/ChConstants.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChSparsityPatternLearner.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

#include <Spectra/KrylovSchurGEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseRegularInverse.h>
#include <Spectra/GenEigsBase.h>

#include <unsupported/Eigen/SparseExtra>  //TODO: remove after debug

using namespace Spectra;
using namespace Eigen;

namespace chrono {

namespace modal {

void BuildGeneralizedEigenProblemMatrices(ChAssembly& assembly,
                                          ChSystemDescriptor& temp_descriptor,
                                          ChSparseMatrix& A,
                                          ChSparseMatrix& B,
                                          int n_vars) {
    // Stiffness matrix
    assembly.LoadKRMMatrices(+1.0, 0.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(A, 0, 0);

    // Mass matrix
    assembly.LoadKRMMatrices(0.0, 0.0, 1.0);
    temp_descriptor.SetMassFactor(1.0);
    temp_descriptor.PasteMassKRMMatrixInto(B, 0, 0);

    // Constraint Jacobian
    assembly.LoadConstraintJacobians();
    temp_descriptor.PasteConstraintsJacobianMatrixInto(A, n_vars, 0);
    temp_descriptor.PasteConstraintsJacobianMatrixTransposedInto(A, 0, n_vars);
}

int ChModalSolverUndamped::Solve(ChAssembly& assembly,
                                 ChMatrixDynamic<double>& eigvects,
                                 ChVectorDynamic<double>& eigvals,
                                 ChVectorDynamic<double>& freq) const {
    m_timer_matrix_assembly.start();

    ChSystemDescriptor sysd;
    ChSystemDescriptor temp_descriptor;

    temp_descriptor.BeginInsertion();
    assembly.InjectVariables(temp_descriptor);
    assembly.InjectKRMMatrices(temp_descriptor);
    assembly.InjectConstraints(temp_descriptor);
    temp_descriptor.EndInsertion();

    temp_descriptor.UpdateCountsAndOffsets();

    // Generate the A and B in state space
    int n_vars = temp_descriptor.CountActiveVariables();
    int n_constr = temp_descriptor.CountActiveConstraints();

    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]

    // B  =  [  M     0  ]
    //       [  0     0  ]

    ChSparsityPatternLearner A_spl(n_vars + n_constr, n_vars + n_constr);
    ChSparsityPatternLearner B_spl(n_vars + n_constr, n_vars + n_constr);
    BuildGeneralizedEigenProblemMatrices(assembly, temp_descriptor, A_spl, B_spl, n_vars);

    ChSparseMatrix A(n_vars + n_constr, n_vars + n_constr);
    ChSparseMatrix B(n_vars + n_constr, n_vars + n_constr);
    A_spl.Apply(A);
    B_spl.Apply(B);

    A.setZeroValues();
    B.setZeroValues();

    // Eigen::saveMarket (const SparseMatrixType &mat, const std::string &filename, int sym=0)

    BuildGeneralizedEigenProblemMatrices(assembly, temp_descriptor, A, B, n_vars);

    // Scale constraints matrix
    double scaling = 1.0;
    if (m_scaleCq) {
        scaling = 0.0;
        for (int k = 0; k < A.outerSize(); ++k) {
            for (ChSparseMatrix::InnerIterator it(A, k); it; ++it) {
                if (it.row() < n_vars && it.col() == it.row()) {
                    scaling += it.valueRef();
                }
            }
        }
        scaling = -scaling / n_vars;
    }

    // A scaling
    for (unsigned int k = 0; k < n_vars + n_constr; ++k) {
        for (ChSparseMatrix::InnerIterator it(A, k); it; ++it) {
            it.valueRef() *= -scaling;
        }
    }

    A.makeCompressed();
    B.makeCompressed();


    std::list<std::pair<int, double>> eig_requests;
    for (int i = 0; i < m_freq_spans.size(); i++) {
        eig_requests.push_back(std::make_pair(m_freq_spans[i].nmodes, m_solver.GetOptimalShift(m_freq_spans[i].freq)));
    }

    m_timer_matrix_assembly.stop();

    m_timer_eigen_solver.start();
    int found_eigs = modal::Solve<>(m_solver, A, B, eigvects, eigvals, eig_requests, m_clip_position_coords ? n_vars : A.rows());
    m_timer_eigen_solver.stop();

    m_timer_solution_postprocessing.start();
    m_solver.GetNaturalFrequencies(eigvals, freq);
    m_timer_solution_postprocessing.stop();

    return found_eigs;
    
}

}  // end namespace modal

}  // end namespace chrono
