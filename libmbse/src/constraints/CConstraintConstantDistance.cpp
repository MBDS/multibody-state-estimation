/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/CConstraintConstantDistance.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

/* -------------------------------------------------------------------
				  CConstraintConstantDistance
   -------------------------------------------------------------------*/
void CConstraintConstantDistance::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless constraint added between two fixed points!");
}

void CConstraintConstantDistance::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	// Update Phi[i]
	// ----------------------------------
	const double dist2 = square(Ax) + square(Ay);
	const double PhiVal = dist2 - square(length);
	arm.Phi_[idx_constr_[0]] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] = 2 * Ax * Adotx + 2 * Ay * Adoty;

	auto& j = jacob.at(0);	// 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	set(j.dPhi_dx[0], -2 * Ax);
	set(j.dPhi_dy[0], -2 * Ay);
	set(j.dPhi_dx[1], 2 * Ax);
	set(j.dPhi_dy[1], 2 * Ay);

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	set(j.dot_dPhi_dx[0], -2 * Adotx);
	set(j.dot_dPhi_dy[0], -2 * Adoty);
	set(j.dot_dPhi_dx[1], 2 * Adotx);
	set(j.dot_dPhi_dy[1], 2 * Adoty);

	// Update Jacobian \{dPhiq*dq}_{\dot{q}}(i,:)
	// --------------------------------------
	set(j.dPhiqdq_dx[0], -2 * Adotx);
	set(j.dPhiqdq_dy[0], -2 * Adoty);
	set(j.dPhiqdq_dx[1], 2 * Adotx);
	set(j.dPhiqdq_dy[1], 2 * Adoty);
}