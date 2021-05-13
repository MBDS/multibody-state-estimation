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

/** \file  This file includes common stuff required by all files in the MBS C++
 * lib.
 */

#include <variant>
#include <memory>  // for auto_ptr
#include <mrpt/poses/CPose3D.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/img/TColor.h>
#include <mrpt/system/CTimeLogger.h>

#include <Eigen/Dense>	// provided by MRPT or standalone
#if EIGEN_VERSION_AT_LEAST(3, 1, 0)
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#else
#error "This library needs Eigen 3.1.0 or newer!"
#endif
// TODO: make this optional?
#include <cholmod.h>
#include <klu.h>

namespace mbse
{
extern mrpt::system::CTimeLogger& timelog();

using dof_index_t = std::size_t;
using point_index_t = std::size_t;

constexpr dof_index_t INVALID_DOF = static_cast<dof_index_t>(-1);
constexpr point_index_t INVALID_POINT_INDEX = static_cast<point_index_t>(-1);

/** Each of the 2D points in a ModelDefinition */
struct Point2
{
	Point2() = default;

	mrpt::math::TPoint2D coords{0, 0};
	/** true: a fixed point, false: generalized coordinates in `q` */
	bool fixed = false;
};

enum class PointDOF : uint8_t
{
	X = 0,
	Y = 1,
	Z = 2
};

/** Return x, y, z */
char dof2letter(const PointDOF p);

/** Degree of Freedom (DOF) info for Natural Coordinates DOFs */
struct NaturalCoordinateDOF
{
	/** Point index in the MBS model */
	point_index_t point_index = INVALID_POINT_INDEX;

	/** DOF, from those of the selected point (0:x, 1:y, 2:z) */
	PointDOF point_dof;

	NaturalCoordinateDOF(point_index_t _point_index, PointDOF _point_dof)
		: point_index(_point_index), point_dof(_point_dof)
	{
	}
};

/**
 * \code
 *        pt2
 *       /
 *      /
 *     / ) angle
 * pt0 ----------- pt1
 * \endcode
 */
struct RelativeAngleDOF
{
	point_index_t point_idx0 = INVALID_POINT_INDEX;
	point_index_t point_idx1 = INVALID_POINT_INDEX;
	point_index_t point_idx2 = INVALID_POINT_INDEX;

	RelativeAngleDOF() = default;
	RelativeAngleDOF(point_index_t i0, point_index_t i1, point_index_t i2)
		: point_idx0(i0), point_idx1(i1), point_idx2(i2)
	{
	}
};

/**
 * \code
 *        pt1
 *       /
 *      /
 *     / ) angle
 * pt0 -----------> +X axis
 * \endcode
 */
struct RelativeAngleAbsoluteDOF
{
	point_index_t point_idx0 = INVALID_POINT_INDEX;
	point_index_t point_idx1 = INVALID_POINT_INDEX;

	RelativeAngleAbsoluteDOF() = default;
	RelativeAngleAbsoluteDOF(point_index_t i0, point_index_t i1)
		: point_idx0(i0), point_idx1(i1)
	{
	}
};

struct RelativeDistanceDOF
{
	point_index_t point_idx0 = INVALID_POINT_INDEX;
	point_index_t point_idx1 = INVALID_POINT_INDEX;

	RelativeDistanceDOF() = default;
	RelativeDistanceDOF(point_index_t i0, point_index_t i1)
		: point_idx0(i0), point_idx1(i1)
	{
	}
};

using RelativeDOF = std::variant<
	std::monostate, RelativeAngleDOF, RelativeAngleAbsoluteDOF,
	RelativeDistanceDOF>;

struct Point2ToDOF
{
	/** The 0-based index in vector q of the `x` and `y` coordinate of this
	 * point (INVALID_DOF if it's a fixed point) */
	dof_index_t dof_x = INVALID_DOF, dof_y = INVALID_DOF;

	Point2ToDOF() = default;
};

struct CompressedRowSparseMatrix
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
		ASSERT_ABOVE_(getNumRows(), 0U);
		ASSERT_ABOVE_(getNumCols(), 0U);
		M.resize(getNumRows(), getNumCols());
		M.fill(0);
		for (size_t row = 0; row < matrix.size(); row++)
			for (const auto& row_val : matrix[row])
			{
				ASSERT_BELOW_(row_val.first, ncols);
				M(row, row_val.first) = row_val.second;
			}
	}

	/** Create a dense version of this sparse matrix */
	template <class MATRIX = Eigen::MatrixXd>
	MATRIX asDense() const
	{
		MATRIX m;
		asDense(m);
		return m;
	}
};

}  // namespace mbse
