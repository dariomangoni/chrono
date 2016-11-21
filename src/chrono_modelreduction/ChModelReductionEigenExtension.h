

//#include "core/ChCSR3Matrix.h"


//void swap(CompressedStorage& other)
//{
//    std::swap(m_values, other.m_values);
//    std::swap(m_indices, other.m_indices);
//    std::swap(m_size, other.m_size);
//    std::swap(m_allocatedSize, other.m_allocatedSize);
//}

template<typename ChCSR3MatrixIMPORTED>
void LoadChCSR3Matrix(const ChCSR3MatrixIMPORTED& ch_mat)
{
    assert(ch_mat.IsRowMajor() == IsRowMajor);
    const_cast<ChCSR3MatrixIMPORTED&>(ch_mat).Compress();
    m_outerSize = IsRowMajor? ch_mat.GetNumRows() : ch_mat.GetNumColumns();     //m_outerSize = leading_dimension;
    m_innerSize = IsRowMajor ? ch_mat.GetNumColumns() : ch_mat.GetNumRows();
    std::free(m_outerIndex);
    m_outerIndex = ch_mat.GetCSR_LeadingIndexArray();

    if (m_innerNonZeros)
    {
        std::free(m_innerNonZeros);
        m_innerNonZeros = 0;
    }

    memset(m_outerIndex, 0, (m_outerSize + 1) * sizeof(Index));

    m_outerIndex    = ch_mat.GetCSR_LeadingIndexArray();
    m_data.index(0) = *ch_mat.GetCSR_TrailingIndexArray();
    m_data.value(0) = *ch_mat.GetCSR_ValueArray();

}


