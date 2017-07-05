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
//
// Interior-Point Utilities
//
// =============================================================================

#ifndef CHIPUTILS_H
#define CHIPUTILS_H
#include "core/ChMatrix.h"
#include "core/ChSparseMatrix.h"
#include <functional>

using namespace chrono;

void PrintMatrix(ChMatrix<>& matrice) {
    for (auto i = 0; i < matrice.GetRows(); i++)
    {
        for (auto j = 0; j < matrice.GetColumns(); j++)
        {
            printf("%.1f ", matrice.GetElement(i, j));
        }
        printf("\n");
    }
}

void PrintCSR3Matrix(ChSparseMatrix& matrice) {
    for (auto i = 0; i < matrice.GetNumRows(); i++)
    {
        for (auto j = 0; j < matrice.GetNumColumns(); j++)
        {
            printf("%3.1f ", matrice.GetElement(i, j));
        }
        printf("\n");
    }
}

template <class matrix>
void ExportArrayToFile(matrix mat, std::string filepath, int precision = 12) {
    std::ofstream ofile;
    ofile.open(filepath);
    ofile << std::scientific << std::setprecision(precision);

    for (auto row_sel = 0; row_sel < mat.GetRows(); row_sel++)
    {
        for (auto col_sel = 0; col_sel < mat.GetColumns(); col_sel++)
        {
            ofile << mat.GetElement(row_sel, col_sel);
        }

        ofile << std::endl;
    }

    ofile.close();
}

void ImportArrayFromFile(ChMatrix<>& output_mat, std::string filename) {
    std::ifstream my_file;
    my_file.open(filename);

    double temp;
    int row_sel = -1;
    for (row_sel = 0; row_sel < output_mat.GetRows(); row_sel++)
    {
        my_file >> temp;
        output_mat.SetElement(row_sel, 0, temp);
    }
    my_file.close();
}




#endif