/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/mbse-common.h>

namespace mbse
{
/** Insert a generic Eigen matrix into a triplet matrix in CHOLMOD format.
 */
template <class MATRIX>
void insert_submatrix_in_triplet(
	cholmod_triplet* tri, const size_t row, const size_t col, const MATRIX& m)
{
	// Insert the entire matrix M or just the upper part?
	if (tri->stype == 0 /* tri is unsymmetric */ || row != col)
	{
		// tri is unsymmetric OR we're out of its diagonal blocks:
		if (tri->stype > 0)
		{
			ASSERTDEB_(col > row);
		}
		else if (tri->stype < 0)
		{
			ASSERTDEB_(row > col);
		}

		// Insert the entire submatrix:
		for (int r = 0; r < m.rows(); r++)
		{
			for (int c = 0; c < m.cols(); c++)
			{
				static_cast<int*>(tri->i)[tri->nnz] = r + row;
				static_cast<int*>(tri->j)[tri->nnz] = c + col;
				static_cast<double*>(tri->x)[tri->nnz] = m.coeff(r, c);
				tri->nnz++;
			}
		}
	}
	else
	{
		// tri is symmetric AND we're on a diagonal block:
		ASSERTDEB_(m.cols() == m.rows());

		if (tri->stype > 0)
		{
			// Insert the upper triangle:
			for (int r = 0; r < m.rows(); r++)
			{
				for (int c = r; c < m.cols(); c++)
				{
					static_cast<int*>(tri->i)[tri->nnz] = r + row;
					static_cast<int*>(tri->j)[tri->nnz] = c + col;
					static_cast<double*>(tri->x)[tri->nnz] = m.coeff(r, c);
					tri->nnz++;
				}
			}
		}
		else
		{
			// Insert the lower triangle:
			for (int r = 0; r < m.rows(); r++)
			{
				for (int c = 0; c <= r; c++)
				{
					static_cast<int*>(tri->i)[tri->nnz] = r + row;
					static_cast<int*>(tri->j)[tri->nnz] = c + col;
					static_cast<double*>(tri->x)[tri->nnz] = m.coeff(r, c);
					tri->nnz++;
				}
			}
		}
	}
}  // end insert_submatrix_in_triplet()

/** Saves a cholmod matrix to a MATLAB file readable with A=load('file');
 * M=spconvert(A); \return false on any error
 */
bool save_matrix(cholmod_sparse* tri, const char* filename, cholmod_common* c);
bool save_matrix_dense(
	cholmod_sparse* tri, const char* filename, cholmod_common* c);

template <typename T, class MATRIX>
void insert_submatrix(
	Eigen::SparseMatrix<T>& A, const size_t row, const size_t col,
	const MATRIX& m)
{
	for (int r = 0; r < m.rows(); r++)
		for (int c = 0; c < m.cols(); c++)
			A.insert(row + r, col + c) = m.coeff(r, c);
}

template <typename T, class MATRIX>
void insert_submatrix(
	std::vector<Eigen::Triplet<T>>& tri, const size_t row, const size_t col,
	const MATRIX& m)
{
	for (int r = 0; r < m.rows(); r++)
		for (int c = 0; c < m.cols(); c++)
			tri.push_back(Eigen::Triplet<T>(row + r, col + c, m.coeff(r, c)));
}

/** Removes columns of the matrix.
 * This "unsafe" version assumes indices sorted in ascending order. */
template <class MATRIX>
void unsafeRemoveColumns(MATRIX& m, const std::vector<std::size_t>& idxs)
{
	std::size_t k = 1;
	const auto nR = m.rows();
	for (auto it = idxs.rbegin(); it != idxs.rend(); ++it, ++k)
	{
		const auto nC = m.cols() - *it - k;
		if (nC > 0)
			m.block(0, *it, nR, nC) = m.block(0, *it + 1, nR, nC).eval();
	}
	m.conservativeResize(nR, m.cols() - idxs.size());
}

/** Removes columns of the matrix. Indices may be unsorted and duplicated */
template <class MATRIX>
void removeColumns(MATRIX& m, const std::vector<std::size_t>& idxsToRemove)
{
	std::vector<std::size_t> idxs = idxsToRemove;
	std::sort(idxs.begin(), idxs.end());
	auto itEnd = std::unique(idxs.begin(), idxs.end());
	idxs.resize(itEnd - idxs.begin());
	unsafeRemoveColumns(m, idxs);
}

/** Extract a subset from a vector (e.g. to extract independent coordinates z
 * from q) */
template <class VECTOR>
VECTOR subset(VECTOR& v, const std::vector<std::size_t>& idxsToExtract)
{
	VECTOR ret;
	ret.resize(idxsToExtract.size());
	for (size_t i = 0; i < idxsToExtract.size(); i++)
		ret(i) = v(idxsToExtract[i]);
	return ret;
}

/** Overwrite part of a vector (e.g. to set independent coordinates z in q) */
template <class SUBVECTOR, class VECTOR>
void overwrite_subset(
	VECTOR& v, const SUBVECTOR& subvector,
	const std::vector<std::size_t>& idxsToOverwrite)
{
	ASSERT_EQUAL_(idxsToOverwrite.size(), subvector.size());
	for (size_t i = 0; i < idxsToOverwrite.size(); i++)
		v(idxsToOverwrite[i]) = subvector(i);
}

}  // namespace mbse
