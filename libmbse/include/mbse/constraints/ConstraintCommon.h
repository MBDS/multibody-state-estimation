/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/AssembledRigidModel.h>
#include <mrpt/core/exceptions.h>
#include <cstdlib>
#include <array>

namespace mbse
{
/** Common data structures and methods for most constraints */
template <
	std::size_t NUM_POINTS, std::size_t NUM_RELATIVE_COORDS = 0,
	std::size_t NUM_JACOB_ROWS = 1>
class ConstraintCommon
{
   public:
	void commonbuildSparseStructures(AssembledRigidModel& arm) const;

	/** Get references to the point coordinates (either fixed or variables in
	 * q) */
	const double& actual_coord(
		const AssembledRigidModel& arm, size_t idx, const PointDOF dof) const
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
		const AssembledRigidModel& arm, size_t idx, const PointDOF dof) const
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

	/** Get references to the point accelerations (either fixed zero or
	 * variables in q) */
	const double& actual_acc(
		const AssembledRigidModel& arm, size_t idx, const PointDOF dof) const
	{
		switch (dof)
		{
			case PointDOF::X:
				return (pointDOFs_[idx].dof_x != INVALID_DOF)
						   ? arm.ddotq_[pointDOFs_[idx].dof_x]
						   : dummy_zero_;
			case PointDOF::Y:
				return (pointDOFs_[idx].dof_y != INVALID_DOF)
						   ? arm.ddotq_[pointDOFs_[idx].dof_y]
						   : dummy_zero_;
			default:
				throw std::invalid_argument("actual_acc(): Invalid dof");
		};
	}

	struct PointRef
	{
		const double &x, &y, &dotx, &doty, &ddotx, &ddoty;

		PointRef(
			const double& X, const double& Y, const double& DOTX,
			const double& DOTY, const double& DDOTX, const double& DDOTY)
			: x(X), y(Y), dotx(DOTX), doty(DOTY), ddotx(DDOTX), ddoty(DDOTY)
		{
		}
		std::string asString() const
		{
			using namespace std::string_literals;
			return "(x,y)=("s + std::to_string(x) + " , "s + std::to_string(y) +
				   ") " + "(vx,vy)=("s + std::to_string(dotx) + " , "s +
				   std::to_string(doty) + ")" + "(ax,ay)=("s +
				   std::to_string(ddotx) + " , "s + std::to_string(ddoty) + ")";
		}
	};

	/** Returns all constant references to X,Y coordinates and velocities of a
	 * point, no matter if the point is fixed (constant) or part of m_q
	 * (generalized coordinates) */
	PointRef actual_coords(const AssembledRigidModel& arm, size_t idx) const
	{
		return {actual_coord(arm, idx, PointDOF::X),
				actual_coord(arm, idx, PointDOF::Y),
				actual_vel(arm, idx, PointDOF::X),
				actual_vel(arm, idx, PointDOF::Y),
				actual_acc(arm, idx, PointDOF::X),
				actual_acc(arm, idx, PointDOF::Y)};
	}

	struct RelCoordRef
	{
		const double &x, &dotx, &ddotx;
		RelCoordRef(const double& X, const double& DOTX, const double& DDOTX)
			: x(X), dotx(DOTX), ddotx(DDOTX)
		{
		}
	};

	/** Gets constant references to a relative coordinate and its velocity.
	 * idxRelativeCoord starts counting at 0 for the first relative coordinate.
	 */
	RelCoordRef actual_rel_coords(
		const AssembledRigidModel& arm, size_t idxRelativeCoord) const
	{
		const auto idx = relativeCoordIndexInQ_.at(idxRelativeCoord);
		return {arm.q_[idx], arm.dotq_[idx], arm.ddotq_[idx]};
	}

   protected:
	ConstraintCommon(
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

	static std::string singlePointDOFAsString(dof_index_t idx)
	{
		if (idx == INVALID_DOF)
			return "fixed";
		else
			return std::to_string(idx);
	}

	std::string pointDOFsAsString() const
	{
		std::string ret;
		for (size_t i = 0; i < NUM_POINTS; i++)
		{
			ret += "point[";
			ret += std::to_string(i);
			ret += "]: (x";
			ret += std::to_string(point_index[i]);
			ret += ",y";
			ret += std::to_string(point_index[i]);
			ret += ") dof_x=";
			ret += singlePointDOFAsString(pointDOFs_[i].dof_x);
			ret += " dof_y=";
			ret += singlePointDOFAsString(pointDOFs_[i].dof_y);
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

		/** Pointers to entries in the sparse Jacobian Phiqq_times_ddq */
		array_dptr_t Phiqq_times_ddq_dx, Phiqq_times_ddq_dy;
		array_relative_dptr_t Phiqq_times_ddq_drel;

		/** Pointers to entries in the sparse Jacobian dotPhiqq_times_dq */
		array_dptr_t dotPhiqq_times_dq_dx, dotPhiqq_times_dq_dy;
		array_relative_dptr_t dotPhiqq_times_dq_drel;

		JacobRowEntries()
		{
			dPhi_dx.fill(nullptr);
			dPhi_dy.fill(nullptr);
			dPhi_drel.fill(nullptr);

			dot_dPhi_dx.fill(nullptr);
			dot_dPhi_dy.fill(nullptr);
			dot_dPhi_drel.fill(nullptr);

			Phiqq_times_ddq_dx.fill(nullptr);
			Phiqq_times_ddq_dy.fill(nullptr);
			Phiqq_times_ddq_drel.fill(nullptr);

			dotPhiqq_times_dq_dx.fill(nullptr);
			dotPhiqq_times_dq_dy.fill(nullptr);
			dotPhiqq_times_dq_drel.fill(nullptr);
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
void ConstraintCommon<NUM_POINTS, NUM_RELATIVE_COORDS, NUM_JACOB_ROWS>::
	commonbuildSparseStructures(AssembledRigidModel& a) const
{
	for (size_t ip = 0; ip < NUM_POINTS; ip++)
	{
		points_[ip] = &a.mechanism_.getPointInfo(point_index[ip]);
		pointDOFs_[ip] = a.getPoints2DOFs()[point_index[ip]];
	}

	// Alloc new rows in the list of constraints:
	for (size_t ic = 0; ic < NUM_JACOB_ROWS; ic++)
	{
		// new Jacobian row index:
		const auto jRow = a.addNewRowToConstraints();
		idx_constr_[ic] = jRow;
		auto& j = jacob[ic];

		auto& Phi_q = a.Phi_q_.matrix;
		auto& dotPhi_q = a.dotPhi_q_.matrix;
		auto& Phiqq_times_ddq = a.Phiqq_times_ddq_.matrix;
		auto& dotPhiqq_times_dq = a.dotPhiqq_times_dq_.matrix;

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

			// Add columns to sparse row in Phiqq_times_ddq:
			j.Phiqq_times_ddq_dx[ip] = &Phiqq_times_ddq[jRow][idxX];
			j.Phiqq_times_ddq_dy[ip] = &Phiqq_times_ddq[jRow][idxY];

			// Add columns to sparse row in dotPhiqq_times_dq:
			j.dotPhiqq_times_dq_dx[ip] = &dotPhiqq_times_dq[jRow][idxX];
			j.dotPhiqq_times_dq_dy[ip] = &dotPhiqq_times_dq[jRow][idxY];
		}

		for (size_t irc = 0; irc < NUM_RELATIVE_COORDS; irc++)
		{
			const auto idxRel = relativeCoordIndexInQ_[irc];

			// Add columns to sparse row in Phi_q:
			j.dPhi_drel[irc] = &Phi_q[jRow][idxRel];

			// Add columns to sparse row in \dot{Phi_q}:
			j.dot_dPhi_drel[irc] = &dotPhi_q[jRow][idxRel];

			// Add columns to sparse row in Phiqq_times_ddq:
			j.Phiqq_times_ddq_drel[irc] = &Phiqq_times_ddq[jRow][idxRel];

			// Add columns to sparse row in dotPhiqq_times_dq:
			j.dotPhiqq_times_dq_drel[irc] = &dotPhiqq_times_dq[jRow][idxRel];
		}
	}
}

}  // namespace mbse
