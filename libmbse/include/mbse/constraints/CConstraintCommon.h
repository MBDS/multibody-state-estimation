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

#include <mbse/CAssembledRigidModel.h>
#include <mrpt/core/exceptions.h>
#include <cstdlib>
#include <array>

namespace mbse
{
/** Common data structures and methods for most constraints */
template <
	std::size_t NUM_POINTS, std::size_t NUM_RELATIVE_COORDS = 0,
	std::size_t NUM_JACOB_ROWS = 1>
class CConstraintCommon
{
   public:
	void commonbuildSparseStructures(CAssembledRigidModel& arm) const;

	/** Get references to the point coordinates (either fixed or variables in
	 * q) */
	const double& actual_coord(
		const CAssembledRigidModel& arm, size_t idx, const PointDOF dof) const
	{
		switch (dof)
		{
			case PointDOF::X:
				return (pointDOFs_[idx].dof_x != INVALID_DOF)
						   ? arm.q_[pointDOFs_[idx].dof_x]
						   : points_[idx]->coords.x;
			case PointDOF::Y:
				return (pointDOFs_[idx].dof_y != INVALID_DOF)
						   ? arm.q_[pointDOFs_[idx].dof_y]
						   : points_[idx]->coords.y;
			default:
				throw std::invalid_argument("actual_coord(): Invalid dof");
		};
	}

	/** Get references to the point velocities (either fixed zero or variables
	 * in q) */
	const double& actual_vel(
		const CAssembledRigidModel& arm, size_t idx, const PointDOF dof) const
	{
		switch (dof)
		{
			case PointDOF::X:
				return (pointDOFs_[idx].dof_x != INVALID_DOF)
						   ? arm.dotq_[pointDOFs_[idx].dof_x]
						   : dummy_zero_;
			case PointDOF::Y:
				return (pointDOFs_[idx].dof_y != INVALID_DOF)
						   ? arm.dotq_[pointDOFs_[idx].dof_y]
						   : dummy_zero_;
			default:
				throw std::invalid_argument("actual_vel(): Invalid dof");
		};
	}

	struct PointRef
	{
		const double &x, &y, &dotx, &doty;

		PointRef(
			const double& X, const double& Y, const double& DOTX,
			const double& DOTY)
			: x(X), y(Y), dotx(DOTX), doty(DOTY)
		{
		}
		std::string asString() const
		{
			using namespace std::string_literals;
			return "(x,y)=("s + std::to_string(x) + " , "s + std::to_string(y) +
				   ") " + "(vx,vy)=("s + std::to_string(dotx) + " , "s +
				   std::to_string(doty) + ")";
		}
	};

	/** Returns all constant references to X,Y coordinates and velocities of a
	 * point, no matter if the point is fixed (constant) or part of m_q
	 * (generalized coordinates) */
	PointRef actual_coords(const CAssembledRigidModel& arm, size_t idx) const
	{
		return {
			actual_coord(arm, idx, PointDOF::X),
			actual_coord(arm, idx, PointDOF::Y),
			actual_vel(arm, idx, PointDOF::X),
			actual_vel(arm, idx, PointDOF::Y)};
	}

	struct RelCoordRef
	{
		const double &x, &dotx;
		RelCoordRef(const double& X, const double& DOTX) : x(X), dotx(DOTX) {}
	};

	/** Gets constant references to a relative coordinate and its velocity.
	 * idxRelativeCoord starts counting at 0 for the first relative coordinate.
	 */
	RelCoordRef actual_rel_coords(
		const CAssembledRigidModel& arm, size_t idxRelativeCoord) const
	{
		return {
			arm.q_[relativeCoordIndexInQ_[idxRelativeCoord]],
			arm.dotq_[relativeCoordIndexInQ_[idxRelativeCoord]]};
	}

   protected:
	CConstraintCommon(
		std::initializer_list<size_t> naturalCoordPointIdxs,
		std::initializer_list<size_t> relativeCoordIndicesInQ =
			std::initializer_list<size_t>())
	{
		points_.fill(nullptr);
		idx_constr_.fill(static_cast<size_t>(-1));
		{
			ASSERT_EQUAL_(naturalCoordPointIdxs.size(), NUM_POINTS);
			size_t idx = 0;
			for (const auto& i : naturalCoordPointIdxs)
				point_index.at(idx++) = i;
		}
		{
			ASSERT_EQUAL_(relativeCoordIndicesInQ.size(), NUM_RELATIVE_COORDS);
			size_t idx = 0;
			for (const auto& i : relativeCoordIndicesInQ)
				relativeCoordIndexInQ_.at(idx++) = i;
		}
	}

	using array_dptr_t = std::array<double*, NUM_POINTS>;
	using array_relative_dptr_t = std::array<double*, NUM_RELATIVE_COORDS>;

	/** Assigned constraint indices, i.e. row indices in the Jacobian dQ_dq
	 */
	mutable std::array<size_t, NUM_JACOB_ROWS> idx_constr_;

	/** Indices of points attached to this constraint: */
	std::array<std::size_t, NUM_POINTS> point_index;

	mutable std::array<const Point2*, NUM_POINTS> points_;

	/** The indices of each point in the state vector "q" */
	mutable std::array<Point2ToDOF, NUM_POINTS> pointDOFs_;

	std::string pointDOFsAsString() const
	{
		std::string ret;
		for (size_t i = 0; i < NUM_POINTS; i++)
		{
			ret += "point[";
			ret += std::to_string(i);
			ret += "]: dof_x=";
			ret += std::to_string(pointDOFs_[i].dof_x);
			ret += " dof_y=";
			ret += std::to_string(pointDOFs_[i].dof_y);
			ret += "\n";
		}
		return ret;
	}

	std::array<size_t, NUM_RELATIVE_COORDS> relativeCoordIndexInQ_;

	struct JacobRowEntries
	{
		/** Pointers to entries in the sparse Jacobian dPhi_dq */
		array_dptr_t dPhi_dx, dPhi_dy;
		array_relative_dptr_t dPhi_drel;  //!< Jacobians wrt rel coords

		/** Pointers to entries in the sparse Jacobian \dot{dPhi_dq} */
		array_dptr_t dot_dPhi_dx, dot_dPhi_dy;
		array_relative_dptr_t dot_dPhi_drel;  //!< Jacobians wrt rel coords

		/** Pointers to entries in the sparse Jacobian d(Phiq*dq)_d{\dot{q}) */
		array_dptr_t dPhiqdq_dx, dPhiqdq_dy;
		array_relative_dptr_t dPhiqdq_drel;	 //!< Jacobians wrt rel coords

		JacobRowEntries()
		{
			dPhi_dx.fill(nullptr);
			dPhi_dy.fill(nullptr);
			dPhi_drel.fill(nullptr);

			dot_dPhi_dx.fill(nullptr);
			dot_dPhi_dy.fill(nullptr);
			dot_dPhi_drel.fill(nullptr);

			dPhiqdq_dx.fill(nullptr);
			dPhiqdq_dy.fill(nullptr);
			dPhiqdq_drel.fill(nullptr);
		}
	};

	const double dummy_zero_ = 0;

	mutable std::array<JacobRowEntries, NUM_JACOB_ROWS> jacob;

	/** Sets a value, if the pointer is !=nullptr */
	static void set(double* ptr, double val)
	{
		if (ptr) *ptr = val;
	}
};

//  ================= Template implementations =============================

template <
	std::size_t NUM_POINTS, std::size_t NUM_RELATIVE_COORDS,
	std::size_t NUM_JACOB_ROWS>
void CConstraintCommon<NUM_POINTS, NUM_RELATIVE_COORDS, NUM_JACOB_ROWS>::
	commonbuildSparseStructures(CAssembledRigidModel& a) const
{
	for (size_t ip = 0; ip < NUM_POINTS; ip++)
	{
		points_[ip] = &a.parent_.getPointInfo(point_index[ip]);
		pointDOFs_[ip] = a.getPoints2DOFs()[point_index[ip]];
	}

	MRPT_TODO("Add Phiqq_times_dq");
	MRPT_TODO("Add d_dotPhiq_ddq_times_dq");

	// Alloc new rows in the list of constraints:
	for (size_t ic = 0; ic < NUM_JACOB_ROWS; ic++)
	{
		// new Jacobian row index:
		const auto jRow = a.addNewRowToConstraints();
		idx_constr_[ic] = jRow;
		auto& j = jacob[ic];

		auto& Phi_q = a.Phi_q_.matrix;
		auto& dotPhi_q = a.dotPhi_q_.matrix;
		auto& dPhiqdq_dq = a.dPhiqdq_dq_.matrix;

		for (size_t ip = 0; ip < NUM_POINTS; ip++)
		{
			// Only for variables, not fixed points
			if (points_[ip]->fixed) continue;

			const auto idxX = pointDOFs_[ip].dof_x;
			const auto idxY = pointDOFs_[ip].dof_y;

			// Add columns to sparse row in Phi_q:
			j.dPhi_dx[ip] = &Phi_q[jRow][idxX];
			j.dPhi_dy[ip] = &Phi_q[jRow][idxY];

			// Add columns to sparse row in \dot{Phi_q}:
			j.dot_dPhi_dx[ip] = &dotPhi_q[jRow][idxX];
			j.dot_dPhi_dy[ip] = &dotPhi_q[jRow][idxY];

			// Add columns to sparse row in d(Phiq*dq)_dq
			j.dPhiqdq_dx[ip] = &dPhiqdq_dq[jRow][idxX];
			j.dPhiqdq_dy[ip] = &dPhiqdq_dq[jRow][idxY];
		}

		for (size_t irc = 0; irc < NUM_RELATIVE_COORDS; irc++)
		{
			const auto idxRel = relativeCoordIndexInQ_[irc];

			// Add columns to sparse row in Phi_q:
			j.dPhi_drel[irc] = &Phi_q[jRow][idxRel];

			// Add columns to sparse row in \dot{Phi_q}:
			j.dot_dPhi_drel[irc] = &dotPhi_q[jRow][idxRel];

			// Add columns to sparse row in d(Phiq*dq)_dq
			j.dPhiqdq_drel[irc] = &dPhiqdq_dq[jRow][idxRel];
		}
	}
}

}  // namespace mbse
