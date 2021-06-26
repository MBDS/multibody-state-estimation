/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/ConstraintRelativeAngle.h>
#include <mrpt/opengl/CSimpleLine.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

void ConstraintRelativeAngle::buildSparseStructures(
	AssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed && points_[2]->fixed),
		"Useless relative coordinate added between three fixed points!");
}

void ConstraintRelativeAngle::update(AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[3] = {
		actual_coords(arm, 0), actual_coords(arm, 1), actual_coords(arm, 2)};

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	// Update Phi[i]
	// ----------------------------------
	MRPT_TODO("Continue!");
	THROW_EXCEPTION("TO DO");
	const double PhiVal = 0;  // XXX
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

	// Update Phiqq_times_ddq
	// ----------------------------------
	MRPT_TODO("Write actual values!");
	set(j.Phiqq_times_ddq_dx[0], 0);
	set(j.Phiqq_times_ddq_dy[0], 0);
	set(j.Phiqq_times_ddq_dx[1], 0);
	set(j.Phiqq_times_ddq_dy[1], 0);

	// Update dotPhiqq_times_dq_dx
	// ----------------------------------
	set(j.dotPhiqq_times_dq_dx[0], 0);
	set(j.dotPhiqq_times_dq_dy[0], 0);
	set(j.dotPhiqq_times_dq_dx[1], 0);
	set(j.dotPhiqq_times_dq_dy[1], 0);
}

void ConstraintRelativeAngle::print(std::ostream& o) const
{
	o << "ConstraintRelativeAngle"
		 "\n";
	o << pointDOFsAsString();
}
