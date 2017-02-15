// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni, Radu Serban
// =============================================================================

#include <algorithm>

#include "chrono/core/ChCSR3Matrix.h"
#include "chrono/core/ChMapMatrix.h"
#include "chrono/solver/ChSystemDescriptor.h"

namespace chrono {

	ChCSR3Matrix::ChCSR3Matrix(int nrows, int ncols, bool row_major_format_on, int nonzeros):
		ChSparseMatrix(nrows,ncols), row_major_format(row_major_format_on)
	{
		// link dimensions to rows and column depending on format
		leading_dimension = row_major_format ? &m_num_rows : &m_num_cols;
		trailing_dimension = row_major_format ? &m_num_cols : &m_num_rows;

		reset_arrays(*leading_dimension, *trailing_dimension, nonzeros);
	}


	/// Put a new element in the matrix following CSR conventions.
	void ChCSR3Matrix::SetElement(int row_sel, int col_sel, double insval, bool overwrite)
	{
		auto lead_sel = row_major_format ? row_sel : col_sel;
		auto trail_sel = row_major_format ? col_sel : row_sel;

		if (insval == 0 && !m_lock)
			return;

		int trail_i;
		for (trail_i = leadIndex[lead_sel]; trail_i<leadIndex[lead_sel+1]; ++trail_i)
		{
			// the requested element DOES NOT exist yet, BUT
			// NO other elements with greater index have been stored yet, SO
			// we can just place the new element here
			if (!initialized_element[trail_i])
			{
				initialized_element[trail_i] = true;
				trailIndex[trail_i] = trail_sel;
				values[trail_i] = insval;
				return;
			}

			// the requested element DOES NOT esist yet AND
			// another element with greater index has already been stored, SO
			// that element has to be pushed further!
			if (trailIndex[trail_i]>trail_sel)
			{
				// insertion needed
                break;
			}

			
			// the requested element already exists
			if (trailIndex[trail_i] == trail_sel)
			{
				(overwrite) ? values[trail_i] = insval : values[trail_i] += insval;
				return;
			}
		}

        // if the loop ends without break means 'row full' (the insertion will also move the other row, for sure!)

		// insertion needed 
		insert(trail_i, lead_sel);
		initialized_element[trail_i] = true;
		trailIndex[trail_i] = trail_sel;
		values[trail_i] = insval;

	}

	double ChCSR3Matrix::GetElement(int row_sel, int col_sel) const
	{
		auto lead_sel = row_major_format ? row_sel : col_sel;
		auto trail_sel = row_major_format ? col_sel : row_sel;

		for (auto trail_i = leadIndex[lead_sel]; trail_i<leadIndex[lead_sel+1]; ++trail_i)
		{
			if (trailIndex[trail_i] == trail_sel)
				return values[trail_i];
		}

		return 0.0;
	}

	double& ChCSR3Matrix::Element(int row_sel, int col_sel)
    {
        auto lead_sel = row_major_format ? row_sel : col_sel;
        auto trail_sel = row_major_format ? col_sel : row_sel;

        int trail_i;
        for (trail_i = leadIndex[lead_sel]; trail_i<leadIndex[lead_sel + 1]; ++trail_i)
        {
            // the requested element DOES NOT exist yet, BUT
            // NO other elements with greater index have been stored yet, SO
            // we can just place the new element here
            if (!initialized_element[trail_i])
            {
                initialized_element[trail_i] = true;
                trailIndex[trail_i] = trail_sel;
                return values[trail_i];
            }

            // the requested element DOES NOT esist yet AND
            // another element with greater index has already been stored, SO
            // that element has to be pushed further!
            if (trailIndex[trail_i]>trail_sel)
            {
                // insertion needed
                break;
            }


            // the requested element already exists
            if (trailIndex[trail_i] == trail_sel)
            {
                return values[trail_i];
            }
        }

        // if the loop ends without break means 'row full' (the insertion will also move the other row, for sure!)

        // insertion needed 
        insert(trail_i, lead_sel);
        initialized_element[trail_i] = true;
        trailIndex[trail_i] = trail_sel;
        return values[trail_i];

    }

	/// The function assures that the matrix will have \a nrows rows and \a ncols columns.
	/// The \a nonzeros_hint value is just a hint!
	/// The matrix can be _fully_ reset, or _partially_ reset (i.e. mantaining the sparsity pattern).
	/// For _partial_ reset the following must apply:
	/// - \a nrows and \a ncols must not differ from the current ones
	/// - \a nonzeros_hint must not be provided (or must equal 0)
	/// - #m_lock must be set
	/// otherwise a _full_ reset will occur.
	void ChCSR3Matrix::Reset(int nrows, int ncols, int nonzeros_hint)
	{
		auto lead_dim_new = row_major_format ? nrows : ncols;
		auto trail_dim_new = row_major_format ? ncols : nrows;


		if (nonzeros_hint == 0 && lead_dim_new == *leading_dimension && trail_dim_new == *trailing_dimension && m_lock)
		{
			std::fill(values.begin(), values.begin() + leadIndex[*leading_dimension] - 1, 0);
		}
		else
		{
			if (nonzeros_hint == 0)
				nonzeros_hint = GetTrailingIndexLength();
			reset_arrays(lead_dim_new, trail_dim_new, nonzeros_hint); // breaks also the sparsity lock
		}
		
		m_num_rows = nrows;
		m_num_cols = ncols;

	}

	int* ChCSR3Matrix::GetCSR_LeadingIndexArray() const
	{
		if (!isCompressed) const_cast<ChCSR3Matrix*>(this)->Compress();
		return const_cast<int*>(leadIndex.data());
	}

	int* ChCSR3Matrix::GetCSR_TrailingIndexArray() const
	{
		if (!isCompressed) const_cast<ChCSR3Matrix*>(this)->Compress();
		return const_cast<int*>(trailIndex.data());
	}

	double* ChCSR3Matrix::GetCSR_ValueArray() const
	{
		if (!isCompressed) const_cast<ChCSR3Matrix*>(this)->Compress();
		return const_cast<double*>(values.data());
	}

	bool ChCSR3Matrix::Compress()
	{
		if (isCompressed)
			return false;

		int trail_i_dest = 0;
		int trail_i = 0;
		for (auto lead_i = 0; lead_i<*leading_dimension; ++lead_i)
		{
			for (; trail_i<leadIndex[lead_i+1]; ++trail_i)
			{
				if (initialized_element[trail_i])
				{
					values[trail_i_dest] = values[trail_i];
					trailIndex[trail_i_dest] = trailIndex[trail_i];
					++trail_i_dest;
				}
			}
			leadIndex[lead_i + 1] = trail_i_dest;
		}

        std::fill(initialized_element.begin(), initialized_element.begin() + leadIndex[*leading_dimension], true);
		isCompressed = true;
		m_lock_broken = false;
		return trail_i_dest!= trail_i;
	}

	void ChCSR3Matrix::Trim()
	{
		trailIndex.resize(leadIndex[*leading_dimension]);
		values.resize(leadIndex[*leading_dimension]);
		initialized_element.resize(leadIndex[*leading_dimension]);

		leadIndex.shrink_to_fit();
		trailIndex.shrink_to_fit();
		values.shrink_to_fit();
		initialized_element.shrink_to_fit();
	}


	void ChCSR3Matrix::Prune(double pruning_threshold)
	{
		int trail_i_dest = 0;
		for (auto lead_i = 0; lead_i<*leading_dimension; ++lead_i)
		{
			for (auto trail_i = leadIndex[lead_i]; trail_i<leadIndex[lead_i+1]; ++trail_i)
			{
				if (initialized_element[trail_i] && std::abs(values[trail_i])>pruning_threshold)
				{
					values[trail_i] = values[trail_i_dest];
					trailIndex[trail_i] = trailIndex[trail_i_dest];
					++trail_i_dest;
				}
			}
			leadIndex[lead_i + 1] = trail_i_dest;
		}
        std::fill(initialized_element.begin(), initialized_element.begin() + leadIndex[*leading_dimension], true);
        m_lock_broken = false;
		isCompressed = true;
	}

	/// Verify if the matrix respects the CSR standard.
	///  3 - warning message: in the row there are no initialized elements
	///  1 - warning message: the matrix is not compressed
	///  0 - all good!
	/// -1 - error message: leadIndex is not strictly ascending
	/// -2 - error message: there's a row that has some an uninitialized element NOT at the end of its space in trailIndex
	/// -4 - error message: trailIndex has not ascending indexes within the rows
	int ChCSR3Matrix::VerifyMatrix() const {
		bool uninitialized_elements_found = false;
		for (int lead_sel = 0; lead_sel < *leading_dimension; lead_sel++) {
			// Check ascending order of leadIndex
			if (leadIndex[lead_sel] >= leadIndex[lead_sel + 1])
			{
				std::cout << "ERROR: leadIndex is not strictly ascending."
					<< " ROW(" << lead_sel << "): " << leadIndex[lead_sel]
					<< "; ROW(" << lead_sel+1 << ") :" << leadIndex[lead_sel + 1]
					<< std::endl;
				return -1;
			}

			bool initialized_elements_found = false;

			int trail_sel = leadIndex[lead_sel + 1];
			while (trail_sel > leadIndex[lead_sel]) {
				trail_sel--;
				if (initialized_element[trail_sel] == false) {
					uninitialized_elements_found = true;
					if (initialized_elements_found)
					{
						std::cout << "ERROR: trailIndex as an invalid not-initialized element at POS: " << trail_sel
							<< "; ROW: " << lead_sel << std::endl;
						return -2;
					}
				}
				else {
					initialized_elements_found = true;
					if (trail_sel > leadIndex[lead_sel] && trailIndex[trail_sel] <= trailIndex[trail_sel - 1])
					{
						std::cout << "ERROR: trailIndex at POS: " << trail_sel-1 <<  " is not greater than the following element."
							<< " COL(" << trail_sel-1 << "): " << trailIndex[trail_sel]
							<< "; COL(" << trail_sel << ") :" << trailIndex[trail_sel]
							<< "; ROW: " << lead_sel << std::endl;
						return -4;
					}
				}
			}
			if (!initialized_elements_found)
			{
				std::cout << "WARNING: no elements in ROW: " << lead_sel << std::endl;
				return 3;
			}
		}

		if (uninitialized_elements_found)
			std::cout << "INFO: the matrix is not compressed" << std::endl;
		else
			std::cout << "OK: matrix verified" << std::endl;

		return (uninitialized_elements_found) ? 1 : 0;
	}

	void ChCSR3Matrix::ImportFromDatFile(std::string path, bool row_major_format_on) {
		std::ifstream a_file, ia_file, ja_file;
		a_file.open(path + "/a.dat");
		ja_file.open(path + "/ja.dat");
		ia_file.open(path + "/ia.dat");

		if (!a_file.is_open())
			assert(0);

		int leadInd_sel;
		for (leadInd_sel = 0; leadInd_sel <= *leading_dimension; leadInd_sel++)
			ia_file >> leadIndex[leadInd_sel];

		Reset(m_num_rows, m_num_cols);

		ia_file.seekg(0);

		for (leadInd_sel = 0; leadInd_sel <= *leading_dimension; leadInd_sel++)
			ia_file >> leadIndex[leadInd_sel];
		leadInd_sel--;

		int trailInd_sel;
		for (trailInd_sel = 0; trailInd_sel < leadIndex[leadInd_sel]; trailInd_sel++) {
			a_file >> values[trailInd_sel];
			ja_file >> trailIndex[trailInd_sel];
		}

		if (trailInd_sel != leadIndex[leadInd_sel])
			assert(0);

		a_file.close();
		ja_file.close();
		ia_file.close();
	}

	void ChCSR3Matrix::ExportToDatFile(std::string filepath, int precision) const {
		std::ofstream a_file, ia_file, ja_file;
		a_file.open(filepath + "/a.dat");
		ja_file.open(filepath + "/ja.dat");
		ia_file.open(filepath + "/ia.dat");
		a_file << std::scientific << std::setprecision(precision);
		ja_file << std::scientific << std::setprecision(precision);
		ia_file << std::scientific << std::setprecision(precision);

		for (int trailInd_sel = 0; trailInd_sel < leadIndex[*leading_dimension]; trailInd_sel++) {
			a_file << values[trailInd_sel] << "\n";
			ja_file << trailIndex[trailInd_sel] << "\n";
		}

		for (int leadInd_sel = 0; leadInd_sel <= *leading_dimension; leadInd_sel++) {
			ia_file << leadIndex[leadInd_sel] << "\n";
		}

		a_file.close();
		ja_file.close();
		ia_file.close();
	}


	/// Acquire information about the sparsity pattern of the elements that _will_ be put into this matrix.
	/// The information is provided by #m_sysd.
	void ChCSR3Matrix::LoadSparsityPattern(ChSparsityPatternLearner& sparsity_learner)
	{
		auto& row_lists = sparsity_learner.GetSparsityPattern();
		*leading_dimension = sparsity_learner.isRowMajor() ? sparsity_learner.GetNumRows() : sparsity_learner.GetNumColumns();
		*trailing_dimension = sparsity_learner.isRowMajor() ? sparsity_learner.GetNumColumns() : sparsity_learner.GetNumRows();
        auto nnz = sparsity_learner.GetNNZ();
        leadIndex.resize(*leading_dimension+1);
		trailIndex.resize(nnz);

		leadIndex[0] = 0;
		for (auto lead_sel = 0; lead_sel<*leading_dimension; ++lead_sel )
		{
			leadIndex[lead_sel + 1] = leadIndex[lead_sel] + static_cast<int>(row_lists[lead_sel].size());
			std::copy(row_lists[lead_sel].begin(), row_lists[lead_sel].end(), trailIndex.begin() + leadIndex[lead_sel]);
		}

		values.resize(nnz);
		initialized_element.assign(nnz, true);
		m_lock_broken = false;
		isCompressed = true;
	}

	void ChCSR3Matrix::distribute_integer_range_on_vector(index_vector_t& vector, int initial_number, int final_number)
	{
		double delta = static_cast<double>(final_number - initial_number) / (vector.size()-1);
		for (auto el_sel = 0; el_sel<vector.size(); el_sel++)
		{
			vector[el_sel] = static_cast<int>(std::round(delta*el_sel));
		}
	}

	void ChCSR3Matrix::reset_arrays(int lead_dim, int trail_dim, int nonzeros)
	{
		
		// break sparsity lock
		m_lock_broken = true;

		// update dimensions (redundant if called from constructor)
		*leading_dimension = lead_dim;
		*trailing_dimension = trail_dim;

		// check if there is at least one element per row
		nonzeros = std::max(nonzeros, lead_dim);

		// allocate the arrays
		leadIndex.resize(*leading_dimension + 1);
        trailIndex.resize(nonzeros);
        values.resize(nonzeros);
        initialized_element.assign(nonzeros, false);

		// make rowIndex span over available space
		distribute_integer_range_on_vector(leadIndex, 0, nonzeros);

		isCompressed = false;

	}

	void ChCSR3Matrix::insert(int& trail_sel, const int& lead_sel)
	{
		isCompressed = false;
		m_lock_broken = true;

		bool OK_also_out_of_row = true; // look for viable positions also in other rows respect to the one selected
		bool OK_also_onelement_rows = false;

		auto trailIndexlength = leadIndex[*leading_dimension];
		int shift_fw = 0; // 0 means no viable position found forward
		int shift_bw = 0; // 0 means no viable position found backward

		//TODO: optimize?
		// look for not-initialized elements FORWARD
		auto lead_sel_fw = lead_sel;
		for (auto trail_i = trail_sel + 1; trail_i < trailIndexlength && (trail_i - trail_sel) < max_shifts; ++trail_i)
		{
			if (!initialized_element[trail_i]) // look for not initialized elements
			{
				// check if it is out of row
				if (!OK_also_out_of_row && trail_i >= leadIndex[lead_sel + 1])
					break;

				// check if it is in 1element row
				if (!OK_also_onelement_rows) // check if we are out of the starting row
				{
					// find the beginning of row that follows the one in which we have found an initialized space
					for (; leadIndex[lead_sel_fw] <= trail_i; ++lead_sel_fw) {}
					if (leadIndex[lead_sel_fw - 1] + 1 >= leadIndex[lead_sel_fw])
						continue;
				}
				shift_fw = trail_i - trail_sel;
				break;
			}
		}

		// look for not-initialized elements BACWARD
		auto lead_sel_bw = lead_sel;
		if (OK_also_out_of_row) // do it only if 'out of row' insertions are allowed
		{
			for (auto trail_i = trail_sel - 1; trail_i >= 0 && (trail_sel - trail_i) < std::min(max_shifts, std::max(shift_fw, 1)); --trail_i)
			{
				if (!initialized_element[trail_i]) // look for not initialized elements
				{
					// check if it is in 1element row
					if (!OK_also_onelement_rows) // check if we are out of the starting row
					{
						// find the beginning of row that follows the one in which we have found an initialized space
						for (; leadIndex[lead_sel_bw] > trail_i; --lead_sel_bw) {}
						if (leadIndex[lead_sel_bw + 1] - 1 <= leadIndex[lead_sel_bw])
							continue;
					}
					shift_bw = trail_i - trail_sel;
					break;
				}
			}
		}
		


		if (shift_bw==0 && shift_fw==0)
		{
			// no viable position found
			// some space has to be make right where trail_sel points
			// meanwhile we give some space also to all the other rows
			// so trail_sel WILL CHANGE

			auto desired_trailIndex_length = static_cast<int>(std::max(GetTrailingIndexLength()*1.2, GetTrailingIndexLength()+1.0));
			auto expansion_factor = 1.5;

			if (desired_trailIndex_length>=trailIndex.size())
			{
				auto new_size = std::max(static_cast<int>(trailIndex.size() * expansion_factor), desired_trailIndex_length);
				index_vector_t trailIndex_new;
				values_vector_t values_new;
				std::vector<bool> initialized_element_new;

				// resize to desired values //TODO: do not initialize elements
                trailIndex_new.resize(new_size);
                values_new.resize(new_size);
                initialized_element_new.assign(new_size, false);

				// copy the array values over new vectors
				copy_and_distribute(trailIndex, values, initialized_element,
									trailIndex_new, values_new, initialized_element_new,
									trail_sel, lead_sel, desired_trailIndex_length - GetTrailingIndexLength());

				// move
				values = std::move(values_new);
				trailIndex = std::move(trailIndex_new);
				initialized_element = std::move(initialized_element_new);
			}
			else
			{
				// resize to desired values
				copy_and_distribute(trailIndex, values, initialized_element,
									trailIndex, values, initialized_element,
									trail_sel, lead_sel, desired_trailIndex_length - GetTrailingIndexLength());
			}

		}
		else
		{
			// shift the elements in order to have a not initialized position where trail_sel points
			// trail_sel WILL CHANGE if backward, WON'T CHANGE if forward
			// WARNING! move_backward is actually forward (towards the end of the array)
			if (shift_bw<0 && -shift_bw<shift_fw)
			{
				if (shift_bw<-1)
				{
					// move is not actually 'moving'; it's just for cleaner code!
					std::move(trailIndex.begin() + trail_sel + shift_bw+1,
							  trailIndex.begin() + trail_sel,
							  trailIndex.begin() + trail_sel + shift_bw);
					std::move(values.begin() + trail_sel + shift_bw + 1,
							  values.begin() + trail_sel,
							  values.begin() + trail_sel + shift_bw);
					std::move(initialized_element.begin() + trail_sel + shift_bw + 1,
							  initialized_element.begin() + trail_sel,
							  initialized_element.begin() + trail_sel + shift_bw);
				}
				

				for (lead_sel_bw = lead_sel; leadIndex[lead_sel_bw] >trail_sel + shift_bw; --lead_sel_bw)
				{
					leadIndex[lead_sel_bw]--;
				}

				trail_sel--;

			}
			else
			{
				// move is not actually 'moving'; it's just for cleaner code!
				std::move_backward(trailIndex.begin() + trail_sel,
								   trailIndex.begin() + trail_sel + shift_fw,
								   trailIndex.begin() + trail_sel + shift_fw+1);
				std::move_backward(values.begin() + trail_sel,
								   values.begin() + trail_sel + shift_fw,
								   values.begin() + trail_sel + shift_fw+1);
				std::move_backward(initialized_element.begin() + trail_sel,
								   initialized_element.begin() + trail_sel + shift_fw,
								   initialized_element.begin() + trail_sel + shift_fw+1);

				for (lead_sel_fw = lead_sel+1; leadIndex[lead_sel_fw] <= trail_sel + shift_fw; ++lead_sel_fw)
				{
					leadIndex[lead_sel_fw]++;
				}
			}

		}
		// let the returning function store the value, takes care of initialize element and so on...

	} // end insert


	void ChCSR3Matrix::copy_and_distribute(const index_vector_t& trailIndex_src,
										   const values_vector_t& values_src,
										   const std::vector<bool>& initialized_element_src,
										   index_vector_t& trailIndex_dest,
										   values_vector_t& values_dest,
										   std::vector<bool>& initialized_element_dest,
										   int& trail_ins, int lead_ins,
										   int storage_augm)
	{
		assert(storage_augm>0);
		assert(trail_ins <= leadIndex[*leading_dimension]);
		assert(leadIndex[*leading_dimension] + storage_augm <= values_dest.capacity());

		double storage_delta = static_cast<double>(storage_augm - 1) / *leading_dimension;

		// start from ending point of source and destination vectors
		auto trail_i_dest = leadIndex[*leading_dimension] + storage_augm;
		auto trail_i_src = leadIndex[*leading_dimension];
		// update leadIndex
		auto lead_i = *leading_dimension;

		while (trail_i_src>0)
		{
			if (trail_i_src == leadIndex[lead_i])
			{
				// evaluate until which element we have to store not-initialized elements
				int fill_up_until = static_cast<int>(std::round(leadIndex[lead_i] + storage_delta*(lead_i - 1) + (((lead_i - 1) >= lead_ins) ? 1 : 0)));

				// insertion point located: make space
				if (trail_i_src == trail_ins && lead_i == lead_ins)
				{
					--trail_i_dest;
					trail_ins = trail_i_dest;
				}

				// update leadIndex
				leadIndex[lead_i] = trail_i_dest;

				--trail_i_dest;
				// fill with not-initialized elements
				while (trail_i_dest >= static_cast<int>(fill_up_until))
				{
					initialized_element_dest[trail_i_dest] = false;
					--trail_i_dest;
				}
				++trail_i_dest;
				--lead_i;

				// insertion point located: make space
				if (trail_i_src == trail_ins && lead_i == lead_ins)
				{
					--trail_i_dest;
					trail_ins = trail_i_dest;
				}

			}
			else
			{
				if (trail_i_src == trail_ins && lead_i == lead_ins)
				{
					--trail_i_dest;
					trail_ins = trail_i_dest;
				}
			}

			// move back and copy elements
			--trail_i_src;
			--trail_i_dest;

			trailIndex_dest[trail_i_dest] = trailIndex_src[trail_i_src];
			initialized_element_dest[trail_i_dest] = initialized_element_src[trail_i_src];
			values_dest[trail_i_dest] = values_src[trail_i_src];
		}
	}

    ChCSR3Matrix& ChCSR3Matrix::operator=(const ChCSR3Matrix& mat_source)
    {
        row_major_format = mat_source.IsRowMajor();
        m_type = mat_source.GetType();
        Reset(mat_source.GetNumRows(), mat_source.GetNumColumns(), mat_source.GetTrailingIndexLength());

        for (auto lead_i = 0; lead_i <= *leading_dimension; ++lead_i)
            leadIndex[lead_i] = mat_source.leadIndex[lead_i];

        for (auto trail_i = 0; trail_i<leadIndex[*leading_dimension]; ++trail_i)
        {
            trailIndex[trail_i] = mat_source.trailIndex[trail_i];
            values[trail_i] = mat_source.values[trail_i];
            initialized_element[trail_i] = mat_source.initialized_element[trail_i];
        }

        isCompressed = mat_source.IsCompressed();
        m_lock_broken = mat_source.m_lock_broken;

        return *this;
    }

    ChCSR3Matrix& ChCSR3Matrix::apply_operator(const ChCSR3Matrix& mat_source, std::function<void(double&, const double&)> f)
    {
        assert(mat_source.GetNumRows() == GetNumRows() && mat_source.GetNumColumns() == GetNumColumns());

        for (auto lead_i = 0; lead_i < *mat_source.leading_dimension; lead_i++)
        {
            for (auto trail_i = mat_source.leadIndex[lead_i]; trail_i < mat_source.trailIndex[lead_i + 1]; trail_i++)
            {
                if (mat_source.initialized_element[trail_i])
                {
                    auto col_i_source = mat_source.IsRowMajor() ? trailIndex[trail_i] : lead_i;
                    auto row_i_source = mat_source.IsRowMajor() ? lead_i : trailIndex[trail_i];

                    f(Element(row_i_source, col_i_source), mat_source.values[trail_i]);

                }
            }
        }

        return *this;
    }

    ChCSR3Matrix& ChCSR3Matrix::operator+=(const ChCSR3Matrix& mat_source)
    {
        return apply_operator(mat_source, [](double& a, const double& b) {a += b; });
    }

    ChCSR3Matrix& ChCSR3Matrix::operator-=(const ChCSR3Matrix& mat_source)
    {
        return apply_operator(mat_source, [](double& a, const double& b) {a -= b; });
    }

    ChCSR3Matrix& ChCSR3Matrix::operator*=(const double coeff)
    {
        assert(coeff != 0 && "ChCSR3Matrix::operator*= cannot accept coeff = zero");

        for (auto trail_i = 0; trail_i < GetTrailingIndexLength(); trail_i++)
            if (initialized_element[trail_i])
                values[trail_i] *= coeff;

        return *this;
    }


    bool ChCSR3Matrix::operator==(const ChCSR3Matrix& mat_source) const
    {
        if (mat_source.GetNumRows() != GetNumRows() || mat_source.GetNumColumns() != GetNumColumns())
            return false;

        for (auto lead_i = 0; lead_i < *leading_dimension; lead_i++)
        {
            auto trail_i_s = mat_source.leadIndex[lead_i];
            auto trail_i_this = this->leadIndex[lead_i];

            while (true)
            {
                if (!initialized_element[trail_i_this] && !mat_source.initialized_element[trail_i_s])
                    break;

                if (trailIndex[trail_i_this] != mat_source.trailIndex[trail_i_s]
                    || values[trail_i_this] != mat_source.values[trail_i_s]
                    || mat_source.initialized_element[trail_i_s] != initialized_element[trail_i_this])
                    return false;
            } //end while

        } //end for

        return true;
    }


    void ChCSR3Matrix::MatrMultiply(const ChMatrix<double>& matB, ChMatrix<double>& mat_out, bool transposeA) const
    {
        transposeA ?
            assert(this->GetNumRows() == matB.GetRows() && "Cannot perform matrix multiplication: wrong dimensions.") :
            assert(this->GetNumColumns() == matB.GetRows() && "Cannot perform matrix multiplication: wrong dimensions.");


        mat_out.Reset(transposeA ? this->GetNumColumns() : this->GetNumRows(), matB.GetColumns());

        for (auto lead_i = 0; lead_i < *leading_dimension; lead_i++)
        {
            for (auto trail_i = leadIndex[lead_i]; trail_i<leadIndex[lead_i + 1] && initialized_element[trail_i]; ++trail_i)
            {
                if (!initialized_element[trail_i])
                    break;

                for (auto col_i = 0; col_i < matB.GetColumns(); col_i++)
                    IsRowMajor() ^ transposeA ? mat_out(lead_i, col_i) += values[trail_i] * matB.GetElement(trailIndex[trail_i], col_i) : mat_out(trailIndex[trail_i], col_i) += values[trail_i] * matB.GetElement(lead_i, col_i);
            }
        }
    }

    
    void ChCSR3Matrix::MatrMultiplyClipped(const ChMatrix<double>& matB, ChMatrix<double>& mat_out,
                                           int start_row_matA, int end_row_matA,
                                           int start_col_matA, int end_col_matA,
                                           int start_row_matB,
                                           int start_row_mat_out,
                                           bool transposeA,
                                           int start_col_matB, int end_col_matB,
                                           int start_col_mat_out) const
    {
        if (end_col_matB == 0)
            end_col_matB = matB.GetColumns() - 1;

        assert(start_row_matA >= 0);
        assert(start_col_matA >= 0);
        assert(start_row_mat_out >= 0);
        assert(start_col_matB >= 0);
        assert(start_row_mat_out >= 0);
        assert(start_col_mat_out >= 0);

        assert(end_col_matB < matB.GetColumns());

        if (!transposeA)
        {
            assert(end_row_matA < GetNumRows());
            assert(end_col_matA < GetNumColumns());
        }
        else
        {
            assert(end_row_matA < GetNumColumns());
            assert(end_col_matA < GetNumRows());
        }

        assert(end_col_matA - start_col_matA + 1 <= matB.GetRows() - start_row_matB && "Not enough rows in matB");
        assert(end_row_matA - start_row_matA + 1 <= mat_out.GetRows() - start_row_mat_out && "Not enough rows in mat_out");
        assert(end_col_matB - start_col_matB + 1 <= mat_out.GetColumns() - start_col_mat_out && "Not enough columns in mat_out");


        auto start_lead_matA = IsRowMajor() ^ transposeA ? start_row_matA : start_col_matA;
        auto end_lead_matA = IsRowMajor() ^ transposeA ? end_row_matA : end_col_matA;
        auto start_trail_matA = IsRowMajor() ^ transposeA ? start_col_matA : start_row_matA;
        auto end_trail_matA = IsRowMajor() ^ transposeA ? end_col_matA : end_row_matA;

        for (auto lead_i = start_lead_matA; lead_i <= end_lead_matA; lead_i++)
        {
            for (auto trail_i = leadIndex[lead_i]; trail_i<leadIndex[lead_i + 1] && initialized_element[trail_i]; ++trail_i)
            {
                if (!initialized_element[trail_i] || trailIndex[trail_i]>end_trail_matA)
                    break;

                if (trailIndex[trail_i] < start_trail_matA)
                    continue;

                for (auto col_i = start_col_matB; col_i <= end_col_matB; col_i++)
                    IsRowMajor() ^ transposeA ?
                    mat_out(start_row_mat_out + lead_i - start_lead_matA,             start_col_mat_out + col_i - start_col_matB) += values[trail_i] * matB.GetElement(start_row_matB + trailIndex[trail_i] - start_col_matA, col_i - start_col_matB) :
                    mat_out(start_row_mat_out + trailIndex[trail_i] - start_col_matA, start_col_mat_out + col_i - start_col_matB) += values[trail_i] * matB.GetElement(start_row_matB + lead_i - start_lead_matA,             col_i - start_col_matB);
            }
        }
    }

    void ChCSR3Matrix::ForEachExistentValueInRange(std::function<void(double*)> func, int start_row, int end_row, int start_col, int end_col)
    {
        auto start_lead = IsRowMajor() ? start_row : start_col;
        auto end_lead = IsRowMajor() ? end_row : end_col;
        auto start_trail = IsRowMajor() ? start_col : start_row;
        auto end_trail = IsRowMajor() ? end_col : end_row;

        for (auto lead_i = start_lead; lead_i <= end_lead; lead_i++)
        {
            for (auto trail_i = leadIndex[lead_i]; trail_i<leadIndex[lead_i + 1] && initialized_element[trail_i]; ++trail_i)
            {
                if (!initialized_element[trail_i] || trailIndex[trail_i]>end_trail)
                    break;

                if (trailIndex[trail_i] < start_trail)
                    continue;

                func(&values[trail_i]);
            }
        }
    }
}  // end namespace chrono

