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
// Authors: Dario Mangoni
// =============================================================================

#ifndef CHCSMATRIXINDEXTYPE_H
#define CHCSMATRIXINDEXTYPE_H

#include "chrono/core/ChCSMatrix.h"
#include "chrono_cvxopt/ChApiCvxopt.h"

#ifdef _DEBUG
#undef _DEBUG
#include <cvxopt.h>
#define _DEBUG
#else
#include <cvxopt.h>
#endif

namespace chrono {
/// \class ChCSMatrixIndexType
/// Class that wraps the ChCSMatrix in order to provide different index types

template <typename index_type = int_t>
class ChCSMatrixIndexType : public chrono::ChCSMatrix {
  public:
    ChCSMatrixIndexType(int nrows, int ncols, bool row_major_format_on, int nonzeros)
        : ChCSMatrix(nrows, ncols, row_major_format_on, nonzeros) {
    }

    virtual ~ChCSMatrixIndexType(){}

    std::vector<index_type> leadIndex_it;
    std::vector<index_type> trailIndex_it;


    bool Compress() override {

        ChCSMatrix::Compress();

        leadIndex_it.resize(m_num_cols + 1);
        trailIndex_it.resize(GetNNZ());

        auto mat_trailIndx = ChCSMatrix::GetCS_TrailingIndexArray();
        auto mat_leadIndx = ChCSMatrix::GetCS_LeadingIndexArray();

        for (auto col_sel = 0; col_sel < m_num_cols + 1; ++col_sel) {
            leadIndex_it[col_sel] = mat_leadIndx[col_sel];
        }

        for (auto nnz_sel = 0; nnz_sel < GetNNZ(); ++nnz_sel) {
            trailIndex_it[nnz_sel] = mat_trailIndx[nnz_sel];
        }

        return true;
    }

    index_type* GetCS_TrailingIndexArray() { return trailIndex_it.data(); }
    index_type* GetCS_LeadingIndexArray() { return leadIndex_it.data(); }

    index_type GetNumRows() const { return static_cast<index_type>(ChCSMatrix::GetNumRows()); }
    index_type GetNumColumns() const { return static_cast<index_type>(ChCSMatrix::GetNumColumns()); }

};

}  // end of namespace chrono

#endif
