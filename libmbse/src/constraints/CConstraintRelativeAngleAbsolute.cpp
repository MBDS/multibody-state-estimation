/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/CConstraintRelativeAngleAbsolute.h>
#include <mrpt/opengl/CSimpleLine.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

void CConstraintRelativeAngleAbsolute::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless relative coordinate added between two fixed points!");
}

void CConstraintRelativeAngleAbsolute::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	MRPT_TODO("Continue!");
	THROW_EXCEPTION("TO DO");

	// Update Phi[i]
	// ----------------------------------
	const double PhiVal = 0;  // XXX
	arm.Phi_[idx_constr_[0]] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] = 0;  // XX

	auto& j = jacob.at(0);	// 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	MRPT_TODO("Continue!");
	if (j.dPhi_dx[0]) *j.dPhi_dx[0] = -2 * Ax;
	if (j.dPhi_dy[0]) *j.dPhi_dy[0] = -2 * Ay;
	if (j.dPhi_dx[1]) *j.dPhi_dx[1] = 2 * Ax;
	if (j.dPhi_dy[1]) *j.dPhi_dy[1] = 2 * Ay;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	MRPT_TODO("Continue!");
	if (j.dot_dPhi_dx[0]) *j.dot_dPhi_dx[0] = 0;
	if (j.dot_dPhi_dy[0]) *j.dot_dPhi_dy[0] = 0;
	if (j.dot_dPhi_dx[1]) *j.dot_dPhi_dx[1] = 0;
	if (j.dot_dPhi_dy[1]) *j.dot_dPhi_dy[1] = 0;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// --------------------------------------
	MRPT_TODO("Continue!");
	if (j.dPhiqdq_dx[0]) *j.dPhiqdq_dx[0] = 0;
	if (j.dPhiqdq_dy[0]) *j.dPhiqdq_dy[0] = 0;
	if (j.dPhiqdq_dx[1]) *j.dPhiqdq_dx[1] = 0;
	if (j.dPhiqdq_dy[1]) *j.dPhiqdq_dy[1] = 0;
}