
#pragma once

/** \file  This file includes common stuff required by all files in the MBS C++
 * lib.
 */

#include <memory>  // for auto_ptr
#include <mrpt/poses/CPose3D.h>

#include <mrpt/version.h>
#if MRPT_VERSION >= 0x199
#include <mrpt/img/TColor.h>
#include <mrpt/system/CTimeLogger.h>
using mrpt::img::TColor;
using mrpt::img::TColorf;
using mrpt::system::CTicTac;
using mrpt::system::CTimeLogger;
using mrpt::system::CTimeLoggerEntry;
#else
#include <mrpt/utils/CTimeLogger.h>
#include <mrpt/utils/TColor.h>
using mrpt::utils::CTicTac;
using mrpt::utils::CTimeLogger;
using mrpt::utils::CTimeLoggerEntry;
using mrpt::utils::TColor;
using mrpt::utils::TColorf;
namespace mrpt
{
using mrpt::math::square;
using mrpt::utils::DEG2RAD;
using mrpt::utils::RAD2DEG;
}  // namespace mrpt
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

extern CTimeLogger timelog;

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
	using row_t = std::map<size_t, double>;

	/** Important: Use deque<> to avoid mem reallocations since we use
	 * *pointers* to elements inside here */
	std::deque<row_t> matrix;

	/** Defines the number of rows. */
	void setRowCount(size_t n) { matrix.resize(n); }

	size_t ncols;  //!< The number of cols in a sparse matrix can be set freely
				   //!< by the user (we don't check for columns out of this
				   //!< limit anyway)

	size_t getNumRows() const { return matrix.size(); }
	size_t getNumCols() const { return ncols; }

	/** Create a dense version of this sparse matrix */
	template <class MATRIX>
	void asDense(MATRIX& M) const
	{
		M.resize(getNumRows(), getNumCols());
		M.fill(0);
		for (size_t i = 0; i < matrix.size(); i++)
			for (const auto& row_val : matrix[i])
				M(i, row_val.first) = row_val.second;
	}
};

}  // namespace sparsembs
