/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CConstraintMobileSlider.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

void CConstraintMobileSlider::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!points_[0]->fixed, "Useless constraint added to a fixed point!");
}

void CConstraintMobileSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& px = (pointDOFs_[0].dof_x != INVALID_DOF)
						   ? arm.q_[pointDOFs_[0].dof_x]
						   : points_[0]->coords.x;
	const double& py = (pointDOFs_[0].dof_y != INVALID_DOF)
						   ? arm.q_[pointDOFs_[0].dof_y]
						   : points_[0]->coords.y;

	const double& pxr0 = (pointDOFs_[1].dof_x != INVALID_DOF)
							 ? arm.q_[pointDOFs_[1].dof_x]
							 : points_[1]->coords.x;
	const double& pyr0 = (pointDOFs_[1].dof_y != INVALID_DOF)
							 ? arm.q_[pointDOFs_[1].dof_y]
							 : points_[1]->coords.y;

	const double& pxr1 = (pointDOFs_[2].dof_x != INVALID_DOF)
							 ? arm.q_[pointDOFs_[2].dof_x]
							 : points_[2]->coords.x;
	const double& pyr1 = (pointDOFs_[2].dof_y != INVALID_DOF)
							 ? arm.q_[pointDOFs_[2].dof_y]
							 : points_[2]->coords.y;

	// Get references to point velocities (fixed=>Zero, variables=>their actual
	// members in dotq_):
	const double dummy_zero = 0;

	const double& vx = (pointDOFs_[0].dof_x != INVALID_DOF)
						   ? arm.dotq_[pointDOFs_[0].dof_x]
						   : dummy_zero;
	const double& vy = (pointDOFs_[0].dof_y != INVALID_DOF)
						   ? arm.dotq_[pointDOFs_[0].dof_y]
						   : dummy_zero;

	const double& vxr0 = (pointDOFs_[1].dof_x != INVALID_DOF)
							 ? arm.dotq_[pointDOFs_[1].dof_x]
							 : dummy_zero;
	const double& vyr0 = (pointDOFs_[1].dof_y != INVALID_DOF)
							 ? arm.dotq_[pointDOFs_[1].dof_y]
							 : dummy_zero;

	const double& vxr1 = (pointDOFs_[2].dof_x != INVALID_DOF)
							 ? arm.dotq_[pointDOFs_[2].dof_x]
							 : dummy_zero;
	const double& vyr1 = (pointDOFs_[2].dof_y != INVALID_DOF)
							 ? arm.dotq_[pointDOFs_[2].dof_y]
							 : dummy_zero;

	// Update Phi[i]
	// ----------------------------------
	arm.Phi_[idx_constr_[0]] =
		(pxr1 - pxr0) * (py - pyr0) - (pyr1 - pyr0) * (px - pxr0);

	// Update dotPhi[i] (partial-Phi[i]_partial-t)
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] =
		(vxr1 - vxr0) * (py - pyr0) + (pxr1 - pxr0) * (vy - vyr0) -
		(vyr1 - vyr0) * (px - pxr0) - (pyr1 - pyr0) * (vx - vxr0);

	auto& j = jacob.at(0);  // 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (j.dPhi_dx[0]) *j.dPhi_dx[0] = pyr0 - pyr1;
	if (j.dPhi_dy[0]) *j.dPhi_dy[0] = pxr1 - pxr0;

	if (j.dPhi_dx[1]) *j.dPhi_dx[1] = -py + pyr1;
	if (j.dPhi_dy[1]) *j.dPhi_dy[1] = px - pxr1;

	if (j.dPhi_dx[2]) *j.dPhi_dx[2] = py - pyr0;
	if (j.dPhi_dy[2]) *j.dPhi_dy[2] = -px + pxr0;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (j.dot_dPhi_dx[0]) *j.dot_dPhi_dx[0] = vyr0 - vyr1;
	if (j.dot_dPhi_dy[0]) *j.dot_dPhi_dy[0] = vxr1 - vxr0;

	if (j.dot_dPhi_dx[1]) *j.dot_dPhi_dx[1] = -vy + vyr1;
	if (j.dot_dPhi_dy[1]) *j.dot_dPhi_dy[1] = vx - vxr1;

	if (j.dot_dPhi_dx[2]) *j.dot_dPhi_dx[2] = vy - vyr0;
	if (j.dot_dPhi_dy[2]) *j.dot_dPhi_dy[2] = -vx + vxr0;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// -------------------------------------
	if (j.dPhiqdq_dx[0]) *j.dPhiqdq_dx[0] = vyr0 - vyr1;
	if (j.dPhiqdq_dy[0]) *j.dPhiqdq_dy[0] = vxr1 - vxr0;
	if (j.dPhiqdq_dx[1]) *j.dPhiqdq_dx[1] = vyr1 - vy;
	if (j.dPhiqdq_dy[1]) *j.dPhiqdq_dy[1] = vx - vxr1;
	if (j.dPhiqdq_dx[2]) *j.dPhiqdq_dx[2] = vy - vyr0;
	if (j.dPhiqdq_dy[2]) *j.dPhiqdq_dy[2] = vxr0 - vx;
}
