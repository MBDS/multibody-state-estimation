/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/ConstraintConstantDistance.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

/* -------------------------------------------------------------------
				  ConstraintConstantDistance
   -------------------------------------------------------------------*/
void ConstraintConstantDistance::buildSparseStructures(
	AssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless constraint added between two fixed points!");
}

void ConstraintConstantDistance::update(AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	const double Addotx = p[1].ddotx - p[0].ddotx;
	const double Addoty = p[1].ddoty - p[0].ddoty;

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
	set(j.dPhi_dx[1], +2 * Ax);
	set(j.dPhi_dy[1], +2 * Ay);

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	set(j.dot_dPhi_dx[0], -2 * Adotx);
	set(j.dot_dPhi_dy[0], -2 * Adoty);
	set(j.dot_dPhi_dx[1], +2 * Adotx);
	set(j.dot_dPhi_dy[1], +2 * Adoty);

	// Update Phiqq_times_ddq
	// ----------------------------------
	set(j.Phiqq_times_ddq_dx[0], -2 * Addotx);
	set(j.Phiqq_times_ddq_dy[0], -2 * Addoty);
	set(j.Phiqq_times_ddq_dx[1], +2 * Addotx);
	set(j.Phiqq_times_ddq_dy[1], +2 * Addoty);

	// Update dotPhiqq_times_dq_dx
	// ----------------------------------
	set(j.dotPhiqq_times_dq_dx[0], 0);
	set(j.dotPhiqq_times_dq_dy[0], 0);
	set(j.dotPhiqq_times_dq_dx[1], 0);
	set(j.dotPhiqq_times_dq_dy[1], 0);
}

void ConstraintConstantDistance::print(std::ostream& o) const
{
	o << "ConstraintConstantDistance, L=" << length << "\n";
	o << pointDOFsAsString();
}
