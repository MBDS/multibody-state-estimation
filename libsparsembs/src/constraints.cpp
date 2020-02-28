#include <sparsembs/constraints.h>
#include <sparsembs/CAssembledRigidModel.h>
#include <mrpt/opengl/CSimpleLine.h>

using namespace sparsembs;
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
	m_points[0] = &arm.m_parent.getPointInfo(this->point_index0);
	m_points[1] = &arm.m_parent.getPointInfo(this->point_index1);

	ASSERTMSG_(
		!(m_points[0]->fixed && m_points[1]->fixed),
		"Useless constraint added between two fixed points!");

	m_pointDOFs[0] = arm.getPoints2DOFs()[this->point_index0];
	m_pointDOFs[1] = arm.getPoints2DOFs()[this->point_index1];

	// Alloc a new row in the list of constraints:
	m_idx_constr = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	if (!m_points[0]->fixed)  // Only for the variables, not fixed points:
	{
		dPhi_dx0 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOFs[0].dof_x];
		dPhi_dy0 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOFs[0].dof_y];
	}
	if (!m_points[1]->fixed)
	{
		dPhi_dx1 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOFs[1].dof_x];
		dPhi_dy1 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOFs[1].dof_y];
	}

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	if (!m_points[0]->fixed)  // Only for the variables, not fixed points:
	{
		dot_dPhi_dx0 =
			&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOFs[0].dof_x];
		dot_dPhi_dy0 =
			&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOFs[0].dof_y];
	}
	if (!m_points[1]->fixed)
	{
		dot_dPhi_dx1 =
			&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOFs[1].dof_x];
		dot_dPhi_dy1 =
			&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOFs[1].dof_y];
	}

	// Add columns to sparse row in d(Phiq*dq)_dq
	// --------------------------------------------
	if (!m_points[0]->fixed)  // Only for the variables, not fixed points:
	{
		dPhiqdq_dx0 =
			&arm.m_dPhiqdq_dq.matrix[m_idx_constr][m_pointDOFs[0].dof_x];
		dPhiqdq_dy0 =
			&arm.m_dPhiqdq_dq.matrix[m_idx_constr][m_pointDOFs[0].dof_y];
	}
	if (!m_points[1]->fixed)
	{
		dPhiqdq_dx1 =
			&arm.m_dPhiqdq_dq.matrix[m_idx_constr][m_pointDOFs[1].dof_x];
		dPhiqdq_dy1 =
			&arm.m_dPhiqdq_dq.matrix[m_idx_constr][m_pointDOFs[1].dof_y];
	}
}

void CConstraintConstantDistance::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& p0x = (m_pointDOFs[0].dof_x != INVALID_DOF)
							? arm.m_q[m_pointDOFs[0].dof_x]
							: m_points[0]->coords.x;
	const double& p0y = (m_pointDOFs[0].dof_y != INVALID_DOF)
							? arm.m_q[m_pointDOFs[0].dof_y]
							: m_points[0]->coords.y;

	const double& p1x = (m_pointDOFs[1].dof_x != INVALID_DOF)
							? arm.m_q[m_pointDOFs[1].dof_x]
							: m_points[1]->coords.x;
	const double& p1y = (m_pointDOFs[1].dof_y != INVALID_DOF)
							? arm.m_q[m_pointDOFs[1].dof_y]
							: m_points[1]->coords.y;

	// Get references to point velocities (fixed=>Zero, variables=>their actual
	// members in m_dotq):
	const double dummy_zero = 0;

	const double& p0dotx = (m_pointDOFs[0].dof_x != INVALID_DOF)
							   ? arm.m_dotq[m_pointDOFs[0].dof_x]
							   : dummy_zero;
	const double& p0doty = (m_pointDOFs[0].dof_y != INVALID_DOF)
							   ? arm.m_dotq[m_pointDOFs[0].dof_y]
							   : dummy_zero;

	const double& p1dotx = (m_pointDOFs[1].dof_x != INVALID_DOF)
							   ? arm.m_dotq[m_pointDOFs[1].dof_x]
							   : dummy_zero;
	const double& p1doty = (m_pointDOFs[1].dof_y != INVALID_DOF)
							   ? arm.m_dotq[m_pointDOFs[1].dof_y]
							   : dummy_zero;

	const double Ax = p1x - p0x;
	const double Ay = p1y - p0y;

	const double Adotx = p1dotx - p0dotx;
	const double Adoty = p1doty - p0doty;

	// Update Phi[i]
	// ----------------------------------
	const double dist2 = square(Ax) + square(Ay);
	const double PhiVal = dist2 - square(length);
	arm.m_Phi[m_idx_constr] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.m_dotPhi[m_idx_constr] = 2 * Ax * Adotx + 2 * Ay * Adoty;

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
	m_Delta = line_pt[1] - line_pt[0];
	ASSERTMSG_(
		m_Delta.norm() > 0,
		"Error: the two constraint points must define a line but they "
		"coincide");

	m_point = &arm.m_parent.getPointInfo(this->point_index);

	ASSERTMSG_(!m_point->fixed, "Useless constraint added to a fixed point!");

	m_pointDOF = arm.getPoints2DOFs()[this->point_index];

	// Alloc a new row in the list of constraints:
	m_idx_constr = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	dPhi_dx0 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOF.dof_x];
	dPhi_dy0 = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOF.dof_y];

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	dot_dPhi_dx0 = &arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOF.dof_x];
	dot_dPhi_dy0 = &arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOF.dof_y];
}

void CConstraintFixedSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& px = arm.m_q[m_pointDOF.dof_x];
	const double& py = arm.m_q[m_pointDOF.dof_y];

	// Get references to point velocities (members in m_dotq):
	const double& pdotx = arm.m_dotq[m_pointDOF.dof_x];
	const double& pdoty = arm.m_dotq[m_pointDOF.dof_y];

	// Update Phi[i] = Ax * (py-y0) - Ay * (px-x0)
	// --------------------------------------------
	const double py_y0 = py - line_pt[0].y;
	const double px_x0 = px - line_pt[0].x;
	arm.m_Phi[m_idx_constr] = m_Delta.x * py_y0 - m_Delta.y * px_x0;

	// Update dotPhi[i]
	// ----------------------------------
	arm.m_dotPhi[m_idx_constr] = m_Delta.x * pdoty - m_Delta.y * pdotx;

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	*dPhi_dx0 = -m_Delta.y;
	*dPhi_dy0 = m_Delta.x;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	*dot_dPhi_dx0 = 0;
	*dot_dPhi_dy0 = 0;

	// Update Jacobian \{dPhiq*dq}_(dq)(i,:)
	// -------------------------------------
	*dPhiqdq_dx0 = -m_Delta.y;
	*dPhiqdq_dy0 = m_Delta.x;
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
	m_points[0] = &arm.m_parent.getPointInfo(this->point_index);
	m_points[1] = &arm.m_parent.getPointInfo(this->ref_pts[0]);
	m_points[2] = &arm.m_parent.getPointInfo(this->ref_pts[1]);

	ASSERTMSG_(
		!m_points[0]->fixed, "Useless constraint added to a fixed point!");

	m_pointDOF[0] = arm.getPoints2DOFs()[this->point_index];
	m_pointDOF[1] = arm.getPoints2DOFs()[this->ref_pts[0]];
	m_pointDOF[2] = arm.getPoints2DOFs()[this->ref_pts[1]];

	// Alloc a new row in the list of constraints:
	m_idx_constr = arm.addNewRowToConstraints();

	// Add columns to sparse row in Phi_q:
	// --------------------------------------------
	for (int i = 0; i < 3; i++)
	{
		if (!m_points[i]->fixed)  // Only for the variables, not fixed points:
		{
			dPhi_dx[i] = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOF[i].dof_x];
			dPhi_dy[i] = &arm.m_Phi_q.matrix[m_idx_constr][m_pointDOF[i].dof_y];
		}
	}

	// Add columns to sparse row in \dot{Phi_q}:
	// --------------------------------------------
	for (int i = 0; i < 3; i++)
	{
		if (!m_points[i]->fixed)  // Only for the variables, not fixed points:
		{
			dot_dPhi_dx[i] =
				&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOF[i].dof_x];
			dot_dPhi_dy[i] =
				&arm.m_dotPhi_q.matrix[m_idx_constr][m_pointDOF[i].dof_y];
		}
	}
}

void CConstraintMobileSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	const double& px = (m_pointDOF[0].dof_x != INVALID_DOF)
						   ? arm.m_q[m_pointDOF[0].dof_x]
						   : m_points[0]->coords.x;
	const double& py = (m_pointDOF[0].dof_y != INVALID_DOF)
						   ? arm.m_q[m_pointDOF[0].dof_y]
						   : m_points[0]->coords.y;

	const double& pxr0 = (m_pointDOF[1].dof_x != INVALID_DOF)
							 ? arm.m_q[m_pointDOF[1].dof_x]
							 : m_points[1]->coords.x;
	const double& pyr0 = (m_pointDOF[1].dof_y != INVALID_DOF)
							 ? arm.m_q[m_pointDOF[1].dof_y]
							 : m_points[1]->coords.y;

	const double& pxr1 = (m_pointDOF[2].dof_x != INVALID_DOF)
							 ? arm.m_q[m_pointDOF[2].dof_x]
							 : m_points[2]->coords.x;
	const double& pyr1 = (m_pointDOF[2].dof_y != INVALID_DOF)
							 ? arm.m_q[m_pointDOF[2].dof_y]
							 : m_points[2]->coords.y;

	// Get references to point velocities (fixed=>Zero, variables=>their actual
	// members in m_dotq):
	const double dummy_zero = 0;

	const double& vx = (m_pointDOF[0].dof_x != INVALID_DOF)
						   ? arm.m_dotq[m_pointDOF[0].dof_x]
						   : dummy_zero;
	const double& vy = (m_pointDOF[0].dof_y != INVALID_DOF)
						   ? arm.m_dotq[m_pointDOF[0].dof_y]
						   : dummy_zero;

	const double& vxr0 = (m_pointDOF[1].dof_x != INVALID_DOF)
							 ? arm.m_dotq[m_pointDOF[1].dof_x]
							 : dummy_zero;
	const double& vyr0 = (m_pointDOF[1].dof_y != INVALID_DOF)
							 ? arm.m_dotq[m_pointDOF[1].dof_y]
							 : dummy_zero;

	const double& vxr1 = (m_pointDOF[2].dof_x != INVALID_DOF)
							 ? arm.m_dotq[m_pointDOF[2].dof_x]
							 : dummy_zero;
	const double& vyr1 = (m_pointDOF[2].dof_y != INVALID_DOF)
							 ? arm.m_dotq[m_pointDOF[2].dof_y]
							 : dummy_zero;

	// Update Phi[i]
	// ----------------------------------
	arm.m_Phi[m_idx_constr] =
		(pxr1 - pxr0) * (py - pyr0) - (pyr1 - pyr0) * (px - pxr0);

	// Update dotPhi[i]
	// ----------------------------------
	arm.m_dotPhi[m_idx_constr] =
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
}
