/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CConstraintConstantDistance.h>

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
	// Get references to the point coordinates (either fixed or variables in q):
	const double& p0x = (pointDOFs_[0].dof_x != INVALID_DOF)
							? arm.q_[pointDOFs_[0].dof_x]
							: points_[0]->coords.x;
	const double& p0y = (pointDOFs_[0].dof_y != INVALID_DOF)
							? arm.q_[pointDOFs_[0].dof_y]
							: points_[0]->coords.y;

	const double& p1x = (pointDOFs_[1].dof_x != INVALID_DOF)
							? arm.q_[pointDOFs_[1].dof_x]
							: points_[1]->coords.x;
	const double& p1y = (pointDOFs_[1].dof_y != INVALID_DOF)
							? arm.q_[pointDOFs_[1].dof_y]
							: points_[1]->coords.y;

	// Get references to point velocities (fixed=>Zero, variables=>their actual
	// members in dotq_):
	const double dummy_zero = 0;

	const double& p0dotx = (pointDOFs_[0].dof_x != INVALID_DOF)
							   ? arm.dotq_[pointDOFs_[0].dof_x]
							   : dummy_zero;
	const double& p0doty = (pointDOFs_[0].dof_y != INVALID_DOF)
							   ? arm.dotq_[pointDOFs_[0].dof_y]
							   : dummy_zero;

	const double& p1dotx = (pointDOFs_[1].dof_x != INVALID_DOF)
							   ? arm.dotq_[pointDOFs_[1].dof_x]
							   : dummy_zero;
	const double& p1doty = (pointDOFs_[1].dof_y != INVALID_DOF)
							   ? arm.dotq_[pointDOFs_[1].dof_y]
							   : dummy_zero;

	const double Ax = p1x - p0x;
	const double Ay = p1y - p0y;

	const double Adotx = p1dotx - p0dotx;
	const double Adoty = p1doty - p0doty;

	// Update Phi[i]
	// ----------------------------------
	const double dist2 = square(Ax) + square(Ay);
	const double PhiVal = dist2 - square(length);
	arm.Phi_[idx_constr_[0]] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] = 2 * Ax * Adotx + 2 * Ay * Adoty;

	auto& j = jacob.at(0);  // 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (j.dPhi_dx[0]) *j.dPhi_dx[0] = -2 * Ax;
	if (j.dPhi_dy[0]) *j.dPhi_dy[0] = -2 * Ay;
	if (j.dPhi_dx[1]) *j.dPhi_dx[1] = 2 * Ax;
	if (j.dPhi_dy[1]) *j.dPhi_dy[1] = 2 * Ay;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (j.dot_dPhi_dx[0]) *j.dot_dPhi_dx[0] = -2 * Adotx;
	if (j.dot_dPhi_dy[0]) *j.dot_dPhi_dy[0] = -2 * Adoty;
	if (j.dot_dPhi_dx[1]) *j.dot_dPhi_dx[1] = 2 * Adotx;
	if (j.dot_dPhi_dy[1]) *j.dot_dPhi_dy[1] = 2 * Adoty;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// --------------------------------------
	if (j.dPhiqdq_dx[0]) *j.dPhiqdq_dx[0] = -2 * Adotx;
	if (j.dPhiqdq_dy[0]) *j.dPhiqdq_dy[0] = -2 * Adoty;
	if (j.dPhiqdq_dx[1]) *j.dPhiqdq_dx[1] = 2 * Adotx;
	if (j.dPhiqdq_dy[1]) *j.dPhiqdq_dy[1] = 2 * Adoty;
}