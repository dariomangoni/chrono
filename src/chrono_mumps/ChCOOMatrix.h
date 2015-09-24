 #ifndef CHCOOMATRIX_H
 #define CHCOOMATRIX_H

#include "ChApiMumps.h"
#include "chrono/core/ChMatrix.h"
#include "chrono/core/ChSparseMatrix.h"
#include <vector>
#include <core/ChMapMatrix.h>

//TODO: template by one-indexed / zero-indexed

namespace chrono{
	class ChApiMumps ChCOOMatrix : public ChMapMatrix
	{

	private:
        std::vector<int> rowIndex;
        std::vector<int> colIndex;
        std::vector<double> values;
        bool one_indexed = false;

	public:
        explicit ChCOOMatrix(int mat_rows = 3, int mat_cols = 3, bool one_indexed = false) :
            ChMapMatrix(mat_rows, mat_cols), one_indexed(one_indexed) {}

		virtual ~ChCOOMatrix(){}

        void SetOneIndexing(bool val) { one_indexed = val; }

        bool Compress() override {
            dynamic_cast<ChMapMatrix*>(this)->ConvertToCOO(rowIndex, colIndex, values, one_indexed);
            return true;
        }

        int* GetRowIndexAddress() const { return const_cast<int*>(rowIndex.data()); }
        int* GetColIndexAddress() const { return const_cast<int*>(colIndex.data()); }
        double* GetValuesAddress() const { return const_cast<double*>(values.data()); }

	};

}

#endif