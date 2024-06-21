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

void BuildUndampedEigenProblemMatrices(ChAssembly& assembly,
                                          ChSystemDescriptor& temp_descriptor,
                                          ChSparseMatrix& A,
                                          ChSparseMatrix& B,
                                          int n_vars) {
    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]
    // WARNING: Cq and Cq' sign is not flipped yet: the user is expected to flip it during scaling (if any)

    // B  =  [  M     0  ]
    //       [  0     0  ]
    //

    // Stiffness matrix
    assembly.LoadKRMMatrices(-1.0, 0.0, 0.0);
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

}  // end namespace modal

}  // end namespace chrono
