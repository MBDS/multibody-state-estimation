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
template <std::size_t NUM_POINTS, std::size_t NUM_JACOB_ROWS = 1>
class CConstraintCommon
{
   public:
	void commonbuildSparseStructures(CAssembledRigidModel& arm) const;

   protected:
	CConstraintCommon(std::initializer_list<size_t> idxs)
	{
		{
			ASSERT_EQUAL_(idxs.size(), NUM_POINTS);
			size_t idx = 0;
			for (const auto& i : idxs) point_index.at(idx++) = i;
		}
		points_.fill(nullptr);
		idx_constr_.fill(static_cast<size_t>(-1));
	}

	using array_dptr_t = std::array<double*, NUM_POINTS>;

	/** Assigned constraint indices, i.e. row indices in the Jacobian dQ_dq */
	mutable std::array<size_t, NUM_JACOB_ROWS> idx_constr_;

	/** Indices of points attached to this constraint: */
	std::array<std::size_t, NUM_POINTS> point_index;

	mutable std::array<const Point2*, NUM_POINTS> points_;

	/** The indices of each point in the state vector "q" */
	mutable std::array<Point2ToDOF, NUM_POINTS> pointDOFs_;

	struct JacobRowEntries
	{
		/** Pointers to entries in the sparse Jacobian dPhi_dq */
		array_dptr_t dPhi_dx, dPhi_dy;

		/** Pointers to entries in the sparse Jacobian \dot{dPhi_dq} */
		array_dptr_t dot_dPhi_dx, dot_dPhi_dy;

		/** Pointers to entries in the sparse Jacobian d(Phiq*dq)_dq */
		array_dptr_t dPhiqdq_dx, dPhiqdq_dy;

		JacobRowEntries()
		{
			dPhi_dx.fill(nullptr);
			dPhi_dy.fill(nullptr);

			dot_dPhi_dx.fill(nullptr);
			dot_dPhi_dy.fill(nullptr);

			dPhiqdq_dx.fill(nullptr);
			dPhiqdq_dy.fill(nullptr);
		}
	};

	mutable std::array<JacobRowEntries, NUM_JACOB_ROWS> jacob;
};

//  =================== Template implementations =============================

template <std::size_t NUM_POINTS, std::size_t NUM_JACOB_ROWS>
void CConstraintCommon<NUM_POINTS, NUM_JACOB_ROWS>::commonbuildSparseStructures(
	CAssembledRigidModel& a) const
{
	for (size_t ip = 0; ip < NUM_POINTS; ip++)
	{
		points_[ip] = &a.parent_.getPointInfo(point_index[ip]);
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
		auto& dPhiqdq_dq = a.dPhiqdq_dq_.matrix;

		for (size_t ip = 0; ip < NUM_POINTS; ip++)
		{
			// Only for variables, not fixed points
			if (points_[ip]->fixed) continue;

			// Add columns to sparse row in Phi_q:
			j.dPhi_dx[ip] = &Phi_q[jRow][pointDOFs_[ip].dof_x];
			j.dPhi_dy[ip] = &Phi_q[jRow][pointDOFs_[ip].dof_y];

			// Add columns to sparse row in \dot{Phi_q}:
			j.dot_dPhi_dx[ip] = &dotPhi_q[jRow][pointDOFs_[ip].dof_x];
			j.dot_dPhi_dy[ip] = &dotPhi_q[jRow][pointDOFs_[ip].dof_y];

			// Add columns to sparse row in d(Phiq*dq)_dq
			j.dPhiqdq_dx[ip] = &dPhiqdq_dq[jRow][pointDOFs_[ip].dof_x];
			j.dPhiqdq_dy[ip] = &dPhiqdq_dq[jRow][pointDOFs_[ip].dof_y];
		}
	}
}

}  // namespace mbse
