
#pragma once

/** \file  This file includes common stuff required by all files in the MBS C++
 * lib.
 */

#include <memory>  // for auto_ptr
#include <mrpt/math/lightweight_geom_data.h>
#include <mrpt/poses/CPose3D.h>

#include <mrpt/version.h>
#if MRPT_VERSION >= 0x199
#include <mrpt/img/TColor.h>
#include <mrpt/system/CTimeLogger.h>
using mrpt::img::TColor;
using mrpt::system::CTimeLogger;
#else
#include <mrpt/utils/CTimeLogger.h>
#include <mrpt/utils/TColor.h>
using mrpt::utils::CTimeLogger;
using mrpt::utils::TColor;
#endif

// Eigen's include must occur AFTER MRPT's headers:
#include <Eigen/Dense>  // provided by MRPT or standalone

#if EIGEN_VERSION_AT_LEAST(3, 1, 0)
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#else
#error "This library needs Eigen 3.1.0 or newer!"
#endif
#include <cholmod.h>
#include <klu.h>

namespace sparsembs
{
using namespace Eigen;
using namespace mrpt::math;

extern mrpt::utils::CTimeLogger timelog;

/** Each of the points in a CModelDefinition */
struct TMBSPoint
{
	TMBSPoint() : coords(0, 0), fixed(false) {}

	TPoint2D coords;
	bool fixed;  //!< true: a fixed point, false: variables
};

/** Degree of Freedom (DOF) info */
struct TDOF
{
	size_t point_index;
	size_t point_dof;  //!< 0:x, 1:y, 2:z

	TDOF(size_t _point_index, size_t _point_dof)
		: point_index(_point_index), point_dof(_point_dof)
	{
	}
};

#define INVALID_DOF static_cast<size_t>(-1)

struct TPoint2DOF
{
	size_t dof_x;  //!< The 0-based index in vector q of the "x" coordinate of
				   //!< this point (-1 if it's a fixed point)
	size_t dof_y;  //!< The 0-based index in vector q of the "y" coordinate of
				   //!< this point (-1 if it's a fixed point)

	TPoint2DOF() : dof_x(INVALID_DOF), dof_y(INVALID_DOF) {}
};

struct TCompressedRowSparseMatrix
{
	typedef std::map<size_t, double> row_t;
	std::deque<row_t> matrix;  // Use deque<> to avodi mem reallocations

	size_t ncols;  //!< The number of cols in a sparse matrix can be set freely
				   //!< by the user (we don't check for columns out of this
				   //!< limit anyway)

	size_t getNumRows() const { return matrix.size(); }
	size_t getNumCols() const { return ncols; }

	/** Create a dense version of this sparse matrix */
	template <class MATRIX>
	void getAsDense(MATRIX& M) const
	{
		M.resize(getNumRows(), getNumCols());
		M.fill(0);
		for (size_t i = 0; i < matrix.size(); i++)
		{
			for (row_t::const_iterator it = matrix[i].begin();
				 it != matrix[i].end(); ++it)
			{
				M(i, it->first) = it->second;
			}
		}
	}
};

}  // namespace sparsembs
