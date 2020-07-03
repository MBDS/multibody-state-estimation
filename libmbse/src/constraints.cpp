/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints.h>
#include <mbse/CAssembledRigidModel.h>
#include <mrpt/opengl/CSimpleLine.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

/** Virtual destructor (required in any virtual base) */
CConstraintBase::~CConstraintBase() {}

/* -------------------------------------------------------------------
				  CConstraintConstantDistance
   -------------------------------------------------------------------*/
void CConstraintConstantDistance::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	points_[0] = &arm.parent_.getPointInfo(this->point_index0);
	points_[1] = &arm.parent_.getPointInfo(this->point_index1);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless constraint added between two fixed points!");

	pointDOFs_[0] = arm.getPoints2DOFs()[this->point_index0];
	pointDOFs_[1] = arm.getPoints2DOFs()[this->point_index1];

	// Alloc a new row in the list of constraints:
	idx_constr_ = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	if (!points_[0]->fixed)  // Only for the variables, not fixed points:
	{
		dPhi_dx0 = &arm.Phi_q_.matrix[idx_constr_][pointDOFs_[0].dof_x];
		dPhi_dy0 = &arm.Phi_q_.matrix[idx_constr_][pointDOFs_[0].dof_y];
	}
	if (!points_[1]->fixed)
	{
		dPhi_dx1 = &arm.Phi_q_.matrix[idx_constr_][pointDOFs_[1].dof_x];
		dPhi_dy1 = &arm.Phi_q_.matrix[idx_constr_][pointDOFs_[1].dof_y];
	}

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	if (!points_[0]->fixed)  // Only for the variables, not fixed points:
	{
		dot_dPhi_dx0 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOFs_[0].dof_x];
		dot_dPhi_dy0 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOFs_[0].dof_y];
	}
	if (!points_[1]->fixed)
	{
		dot_dPhi_dx1 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOFs_[1].dof_x];
		dot_dPhi_dy1 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOFs_[1].dof_y];
	}

	// Add columns to sparse row in d(Phiq*dq)_dq
	// --------------------------------------------
	if (!points_[0]->fixed)  // Only for the variables, not fixed points:
	{
		dPhiqdq_dx0 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOFs_[0].dof_x];
		dPhiqdq_dy0 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOFs_[0].dof_y];
	}
	if (!points_[1]->fixed)
	{
		dPhiqdq_dx1 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOFs_[1].dof_x];
		dPhiqdq_dy1 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOFs_[1].dof_y];
	}
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
	arm.Phi_[idx_constr_] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_] = 2 * Ax * Adotx + 2 * Ay * Adoty;

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (dPhi_dx0) *dPhi_dx0 = -2 * Ax;
	if (dPhi_dy0) *dPhi_dy0 = -2 * Ay;
	if (dPhi_dx1) *dPhi_dx1 = 2 * Ax;
	if (dPhi_dy1) *dPhi_dy1 = 2 * Ay;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (dot_dPhi_dx0) *dot_dPhi_dx0 = -2 * Adotx;
	if (dot_dPhi_dy0) *dot_dPhi_dy0 = -2 * Adoty;
	if (dot_dPhi_dx1) *dot_dPhi_dx1 = 2 * Adotx;
	if (dot_dPhi_dy1) *dot_dPhi_dy1 = 2 * Adoty;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// --------------------------------------
	if (dPhiqdq_dx0) *dPhiqdq_dx0 = -2 * Adotx;
	if (dPhiqdq_dy0) *dPhiqdq_dy0 = -2 * Adoty;
	if (dPhiqdq_dx1) *dPhiqdq_dx1 = 2 * Adotx;
	if (dPhiqdq_dy1) *dPhiqdq_dy1 = 2 * Adoty;
}

/* -------------------------------------------------------------------
				  CConstraintFixedSlider
   -------------------------------------------------------------------*/
void CConstraintFixedSlider::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	Delta_ = line_pt[1] - line_pt[0];
	ASSERTMSG_(
		Delta_.norm() > 0,
		"Error: the two constraint points must define a line but they "
		"coincide");

	point_ = &arm.parent_.getPointInfo(this->point_index);

	ASSERTMSG_(!point_->fixed, "Useless constraint added to a fixed point!");

	pointDOF_ = arm.getPoints2DOFs()[this->point_index];

	// Alloc a new row in the list of constraints:
	idx_constr_ = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	dPhi_dx0 = &arm.Phi_q_.matrix[idx_constr_][pointDOF_.dof_x];
	dPhi_dy0 = &arm.Phi_q_.matrix[idx_constr_][pointDOF_.dof_y];

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	dot_dPhi_dx0 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOF_.dof_x];
	dot_dPhi_dy0 = &arm.dotPhi_q_.matrix[idx_constr_][pointDOF_.dof_y];

	// Add columns to sparse row in d(Phiq*dq)_dq
	// --------------------------------------------
	dPhiqdq_dx0 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOF_.dof_x];
	dPhiqdq_dy0 = &arm.dPhiqdq_dq_.matrix[idx_constr_][pointDOF_.dof_y];
}

void CConstraintFixedSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& px = arm.q_[pointDOF_.dof_x];
	const double& py = arm.q_[pointDOF_.dof_y];

	// Get references to point velocities (members in dotq_):
	const double& pdotx = arm.dotq_[pointDOF_.dof_x];
	const double& pdoty = arm.dotq_[pointDOF_.dof_y];

	// Update Phi[i] = Ax * (py-y0) - Ay * (px-x0)
	// --------------------------------------------
	const double py_y0 = py - line_pt[0].y;
	const double px_x0 = px - line_pt[0].x;
	arm.Phi_[idx_constr_] = Delta_.x * py_y0 - Delta_.y * px_x0;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_] = Delta_.x * pdoty - Delta_.y * pdotx;

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (dPhi_dx0) *dPhi_dx0 = -Delta_.y;
	if (dPhi_dy0) *dPhi_dy0 = Delta_.x;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (dot_dPhi_dx0) *dot_dPhi_dx0 = 0;
	if (dot_dPhi_dy0) *dot_dPhi_dy0 = 0;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// -------------------------------------
	if (dPhiqdq_dx0) *dPhiqdq_dx0 = -Delta_.y;
	if (dPhiqdq_dy0) *dPhiqdq_dy0 = Delta_.x;
}

/** Creates a 3D representation of the constraint, if applicable (e.g. the line
 * of a fixed slider) \return false if the constraint has no 3D representation
 */
bool CConstraintFixedSlider::get3DRepresentation(
	mrpt::opengl::CRenderizable::Ptr& inout_obj) const
{
	inout_obj = mrpt::opengl::CSimpleLine::Create(
		line_pt[0].x, line_pt[0].y, 0, line_pt[1].x, line_pt[1].y, 0,
		1 /*line width*/);

	return true;
}

/* -------------------------------------------------------------------
				  CConstraintMobileSlider
   -------------------------------------------------------------------*/
void CConstraintMobileSlider::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	points_[0] = &arm.parent_.getPointInfo(this->point_index);
	points_[1] = &arm.parent_.getPointInfo(this->ref_pts[0]);
	points_[2] = &arm.parent_.getPointInfo(this->ref_pts[1]);

	ASSERTMSG_(
		!points_[0]->fixed, "Useless constraint added to a fixed point!");

	pointDOF_[0] = arm.getPoints2DOFs()[this->point_index];
	pointDOF_[1] = arm.getPoints2DOFs()[this->ref_pts[0]];
	pointDOF_[2] = arm.getPoints2DOFs()[this->ref_pts[1]];

	// Alloc a new row in the list of constraints:
	idx_constr_ = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	for (int i = 0; i < 3; i++)
	{
		if (!points_[i]->fixed)  // Only for the variables, not fixed points:
		{
			dPhi_dx[i] = &arm.Phi_q_.matrix[idx_constr_][pointDOF_[i].dof_x];
			dPhi_dy[i] = &arm.Phi_q_.matrix[idx_constr_][pointDOF_[i].dof_y];
		}
	}

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	for (int i = 0; i < 3; i++)
	{
		if (!points_[i]->fixed)  // Only for the variables, not fixed points:
		{
			dot_dPhi_dx[i] =
				&arm.dotPhi_q_.matrix[idx_constr_][pointDOF_[i].dof_x];
			dot_dPhi_dy[i] =
				&arm.dotPhi_q_.matrix[idx_constr_][pointDOF_[i].dof_y];
		}
	}
}

void CConstraintMobileSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& px = (pointDOF_[0].dof_x != INVALID_DOF)
						   ? arm.q_[pointDOF_[0].dof_x]
						   : points_[0]->coords.x;
	const double& py = (pointDOF_[0].dof_y != INVALID_DOF)
						   ? arm.q_[pointDOF_[0].dof_y]
						   : points_[0]->coords.y;

	const double& pxr0 = (pointDOF_[1].dof_x != INVALID_DOF)
							 ? arm.q_[pointDOF_[1].dof_x]
							 : points_[1]->coords.x;
	const double& pyr0 = (pointDOF_[1].dof_y != INVALID_DOF)
							 ? arm.q_[pointDOF_[1].dof_y]
							 : points_[1]->coords.y;

	const double& pxr1 = (pointDOF_[2].dof_x != INVALID_DOF)
							 ? arm.q_[pointDOF_[2].dof_x]
							 : points_[2]->coords.x;
	const double& pyr1 = (pointDOF_[2].dof_y != INVALID_DOF)
							 ? arm.q_[pointDOF_[2].dof_y]
							 : points_[2]->coords.y;

	// Get references to point velocities (fixed=>Zero, variables=>their actual
	// members in dotq_):
	const double dummy_zero = 0;

	const double& vx = (pointDOF_[0].dof_x != INVALID_DOF)
						   ? arm.dotq_[pointDOF_[0].dof_x]
						   : dummy_zero;
	const double& vy = (pointDOF_[0].dof_y != INVALID_DOF)
						   ? arm.dotq_[pointDOF_[0].dof_y]
						   : dummy_zero;

	const double& vxr0 = (pointDOF_[1].dof_x != INVALID_DOF)
							 ? arm.dotq_[pointDOF_[1].dof_x]
							 : dummy_zero;
	const double& vyr0 = (pointDOF_[1].dof_y != INVALID_DOF)
							 ? arm.dotq_[pointDOF_[1].dof_y]
							 : dummy_zero;

	const double& vxr1 = (pointDOF_[2].dof_x != INVALID_DOF)
							 ? arm.dotq_[pointDOF_[2].dof_x]
							 : dummy_zero;
	const double& vyr1 = (pointDOF_[2].dof_y != INVALID_DOF)
							 ? arm.dotq_[pointDOF_[2].dof_y]
							 : dummy_zero;

	// Update Phi[i]
	// ----------------------------------
	arm.Phi_[idx_constr_] =
		(pxr1 - pxr0) * (py - pyr0) - (pyr1 - pyr0) * (px - pxr0);

	// Update dotPhi[i] (partial-Phi[i]_partial-t)
	// ----------------------------------
	arm.dotPhi_[idx_constr_] =
		(vxr1 - vxr0) * (py - pyr0) + (pxr1 - pxr0) * (vy - vyr0) -
		(vyr1 - vyr0) * (px - pxr0) - (pyr1 - pyr0) * (vx - vxr0);

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (dPhi_dx[0]) *dPhi_dx[0] = pyr0 - pyr1;
	if (dPhi_dy[0]) *dPhi_dy[0] = pxr1 - pxr0;

	if (dPhi_dx[1]) *dPhi_dx[1] = -py + pyr1;
	if (dPhi_dy[1]) *dPhi_dy[1] = px - pxr1;

	if (dPhi_dx[2]) *dPhi_dx[2] = py - pyr0;
	if (dPhi_dy[2]) *dPhi_dy[2] = -px + pxr0;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (dot_dPhi_dx[0]) *dot_dPhi_dx[0] = vyr0 - vyr1;
	if (dot_dPhi_dy[0]) *dot_dPhi_dy[0] = vxr1 - vxr0;

	if (dot_dPhi_dx[1]) *dot_dPhi_dx[1] = -vy + vyr1;
	if (dot_dPhi_dy[1]) *dot_dPhi_dy[1] = vx - vxr1;

	if (dot_dPhi_dx[2]) *dot_dPhi_dx[2] = vy - vyr0;
	if (dot_dPhi_dy[2]) *dot_dPhi_dy[2] = -vx + vxr0;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// -------------------------------------
	if (dPhiqdq_dx) *dPhiqdq_dx = vyr0 - vyr1;
	if (dPhiqdq_dy) *dPhiqdq_dy = vxr1 - vxr0;
	if (dPhiqdq_dx0) *dPhiqdq_dx0 = vyr1 - vy;
	if (dPhiqdq_dy0) *dPhiqdq_dy0 = vx - vxr1;
	if (dPhiqdq_dx1) *dPhiqdq_dx1 = vy - vyr0;
	if (dPhiqdq_dy1) *dPhiqdq_dy1 = vxr0 - vx;
}
