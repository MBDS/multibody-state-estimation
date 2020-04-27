/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CAssembledRigidModel.h>
#include <mrpt/opengl.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

const double DEFAULT_GRAVITY[3] = {0, -9.81, 0};

/** Constructor */
CAssembledRigidModel::CAssembledRigidModel(const TSymbolicAssembledModel& armi)
	: m_parent(armi.model)
{
	for (int i = 0; i < 3; i++) m_gravity[i] = DEFAULT_GRAVITY[i];

	const size_t nDOFs = armi.DOFs.size();
	ASSERTMSG_(nDOFs > 0, "Trying to assemble model with 0 DOFs!?");

	m_q.resize(nDOFs);
	m_q.setConstant(0);
	m_dotq.resize(nDOFs);
	m_dotq.setConstant(0);
	m_ddotq.resize(nDOFs);
	m_ddotq.setConstant(0);
	m_DOFs = armi.DOFs;
	m_points2DOFs.resize(armi.model.getPointCount());

	for (size_t i = 0; i < nDOFs; i++)
	{
		const size_t pt_idx = m_DOFs[i].point_index;
		const TMBSPoint& pt = m_parent.getPointInfo(pt_idx);
		switch (m_DOFs[i].point_dof)
		{
			case 0:  // x
				m_q[i] = pt.coords.x;
				m_points2DOFs[pt_idx].dof_x = i;
				break;
			case 1:  // y
				m_q[i] = pt.coords.y;
				m_points2DOFs[pt_idx].dof_y = i;
				break;
			// case 2: //z
			//	m_q[i] = pt.coords.z;
			//	break;
			default:
				THROW_EXCEPTION("Unexpected value for m_DOFs[i].point_dof");
		};
	}

	// Save the number of DOFs as the number of columsn in sparse Jacobians:
	m_Phi_q.ncols = nDOFs;
	m_dotPhi_q.ncols = nDOFs;
	m_dPhiqdq_dq.ncols = nDOFs;

	// Generate constraint equations & create structure of sparse Jacobians:
	// ---------------------------------------------------------------------------------------------------
	const std::vector<CConstraintBase::Ptr>& parent_constraints =
		m_parent.getConstraints();
	const size_t nConst = parent_constraints.size();

	m_constraints.resize(nConst);

	for (size_t i = 0; i < nConst; i++)
	{
		m_constraints[i] = parent_constraints[i];  // Copy smart pointer, and
		// clone object:
		m_constraints[i] = CConstraintBase::Ptr(m_constraints[i]->clone());

		const CConstraintBase* c = m_constraints[i].get();
		ASSERT_(c != NULL);
		c->buildSparseStructures(*this);
	}
}

/** Returns the current gravity aceleration vector, used for the bodies weights
 * (default: [0 -9.81 0]) */
void CAssembledRigidModel::getGravityVector(
	double& gx, double& gy, double& gz) const
{
	gx = m_gravity[0];
	gy = m_gravity[1];
	gz = m_gravity[2];
}

/** Changes the gravity aceleration vector, used for the bodies weights
 * (default: [0 -9.81 0]) */
void CAssembledRigidModel::setGravityVector(
	const double gx, const double gy, const double gz)
{
	m_gravity[0] = gx;
	m_gravity[1] = gy;
	m_gravity[2] = gz;
}

/** Call all constraint objects and command them to update their corresponding
 * parts in the sparse Jacobians */
void CAssembledRigidModel::update_numeric_Phi_and_Jacobians()
{
	// Update numeric values of the constraint Jacobians:
	for (size_t i = 0; i < m_constraints.size(); i++)
		m_constraints[i]->update(*this);
}

/** Returns a 3D visualization of the model */
void CAssembledRigidModel::getAs3DRepresentation(
	mrpt::opengl::CSetOfObjects::Ptr& outObj,
	const CBody::TRenderParams& rp) const
{
	if (!outObj)
		outObj = mrpt::opengl::CSetOfObjects::Create();
	else
		outObj->clear();

	// Render constraints:
	const std::vector<CConstraintBase::Ptr>& parent_constr =
		m_parent.getConstraints();

	const size_t nConstr = parent_constr.size();
	for (size_t i = 0; i < nConstr; i++)
	{
		const CConstraintBase* constr_ptr = parent_constr[i].get();
		mrpt::opengl::CRenderizable::Ptr gl_obj;
		if (constr_ptr->get3DRepresentation(gl_obj))
		{
			// Insert in 3D scene:
			outObj->insert(gl_obj);
		}
	}

	// Render "ground" points:
	if (rp.show_grounds)
	{
		const size_t nPts = m_parent.getPointCount();
		for (size_t i = 0; i < nPts; i++)
		{
			const TMBSPoint& pt = m_parent.getPointInfo(i);
			if (!pt.fixed) continue;
			// This is a fixed point:
			outObj->insert(this->internal_render_ground_point(pt, rp));
		}
	}

	// Render bodies:
	const std::vector<CBody>& parent_bodies = m_parent.getBodies();
	const size_t nBodies = parent_bodies.size();
	m_gl_objects.resize(nBodies);

	for (size_t i = 0; i < nBodies; i++)
	{
		const CBody& b = parent_bodies[i];

		mrpt::opengl::CRenderizable::Ptr gl_obj = b.get3DRepresentation();

		// Insert in 3D scene:
		outObj->insert(gl_obj);

		// And save reference for quickly update the 3D pose in the future
		// during animations:
		m_gl_objects[i] = gl_obj;
	}

	// Place each body in its current pose:
	this->update3DRepresentation(rp);
}

/** Animates a 3D representation of the MBS, previously built in
 * getAs3DRepresentation() \sa getAs3DRepresentation
 */
void CAssembledRigidModel::update3DRepresentation(
	const CBody::TRenderParams& rp) const
{
	const std::vector<CBody>& parent_bodies = m_parent.getBodies();

	const size_t nBodies = parent_bodies.size();
	if (m_gl_objects.size() != nBodies)
	{
		std::cerr << "[CAssembledRigidModel::update3DRepresentation] Warning: "
					 "Opengl model is not initialized.\n";
		return;
	}

	for (size_t i = 0; i < nBodies; i++)
	{
		mrpt::opengl::CRenderizable::Ptr& obj = m_gl_objects[i];
		ASSERT_(obj);

		// Recover the 2D pose from 2 points:
		const CBody& b = parent_bodies[i];
		const size_t i0 = b.points[0];
		const size_t i1 = b.points[1];

		const TMBSPoint& pnt0 = m_parent.getPointInfo(i0);
		const TMBSPoint& pnt1 = m_parent.getPointInfo(i1);

		const TPoint2DOF dof0 = m_points2DOFs[i0];
		const TPoint2DOF dof1 = m_points2DOFs[i1];

		const double& p0x =
			(dof0.dof_x != INVALID_DOF) ? m_q[dof0.dof_x] : pnt0.coords.x;
		const double& p0y =
			(dof0.dof_y != INVALID_DOF) ? m_q[dof0.dof_y] : pnt0.coords.y;

		const double& p1x =
			(dof1.dof_x != INVALID_DOF) ? m_q[dof1.dof_x] : pnt1.coords.x;
		const double& p1y =
			(dof1.dof_y != INVALID_DOF) ? m_q[dof1.dof_y] : pnt1.coords.y;

		const double theta = atan2(p1y - p0y, p1x - p0x);

		obj->setPose(
			mrpt::poses::CPose3D(p0x, p0y, 0, theta, DEG2RAD(0), DEG2RAD(0)));

		// Update transparency:
		if (rp.render_style == CBody::reLine)
		{
			auto gl_line =
				mrpt::ptr_cast<mrpt::opengl::CSetOfObjects>::from(obj)
					->getByClass<mrpt::opengl::CSimpleLine>();
			if (gl_line) gl_line->setColorA_u8(rp.line_alpha);
		}
	}
}

size_t CAssembledRigidModel::addNewRowToConstraints()
{
	const size_t idx = m_Phi.size();
	const size_t m = idx + 1;  // new size

	// Add rows:
	m_Phi.resize(m);
	m_dotPhi.resize(m);

	// Jacobians:
	m_Phi_q.setRowCount(m);
	m_dotPhi_q.setRowCount(m);
	m_dPhiqdq_dq.setRowCount(m);

	return idx;
}

/** Only to be called between objects created from the same symbolic model, this
 * method replicates the state of "o" into "this". */
void CAssembledRigidModel::copyStateFrom(const CAssembledRigidModel& o)
{
#ifdef _DEBUG
	// Security checks, just in case...
	ASSERT_(this->m_q.size() == o.m_q.size());
	ASSERT_(this->m_dotq.size() == o.m_dotq.size());
	ASSERT_(this->m_DOFs.size() == o.m_DOFs.size());
	ASSERT_(this->m_Phi.size() == o.m_Phi.size());
	const double* ptr_q0 = &m_q[0];
	const double* ptr_dotq0 = &m_dotq[0];
#endif

	this->m_q = o.m_q;
	this->m_dotq = o.m_dotq;

#ifdef _DEBUG
	ASSERT_(
		ptr_q0 == &m_q[0]);  // make sure the vectors didn't suffer mem
							 // reallocation, since we save pointers to these!
	ASSERT_(
		ptr_dotq0 ==
		&m_dotq[0]);  // make sure the vectors didn't suffer mem reallocation,
					  // since we save pointers to these!
#endif
}

/** Copies the opengl object from another instance */
void CAssembledRigidModel::copyOpenGLRepresentationFrom(
	const CAssembledRigidModel& o)
{
	this->m_gl_objects = o.m_gl_objects;
}

/** Retrieves the current coordinates of a point, which may include either fixed
 * or variable components */
void CAssembledRigidModel::getPointCurrentCoords(
	const size_t pt_idx, mrpt::math::TPoint2D& pt) const
{
	const TMBSPoint& pt_info = m_parent.getPointInfo(pt_idx);
	const TPoint2DOF& pt_dofs = m_points2DOFs[pt_idx];

	pt.x =
		(pt_dofs.dof_x != INVALID_DOF) ? m_q[pt_dofs.dof_x] : pt_info.coords.x;
	pt.y =
		(pt_dofs.dof_y != INVALID_DOF) ? m_q[pt_dofs.dof_y] : pt_info.coords.y;
}

/** Retrieves the current velocity of a point, which may include either fixed or
 * variable components */
void CAssembledRigidModel::getPointCurrentVelocity(
	const size_t pt_idx, mrpt::math::TPoint2D& vel) const
{
	const TPoint2DOF& pt_dofs = m_points2DOFs[pt_idx];

	vel.x = (pt_dofs.dof_x != INVALID_DOF) ? m_dotq[pt_dofs.dof_x] : 0;
	vel.y = (pt_dofs.dof_y != INVALID_DOF) ? m_dotq[pt_dofs.dof_y] : 0;
}

/** Computes the current coordinates of a point fixed to a given body, given its
 * relative coordinates wrt to system X:pt0->pt1, Y: orthogonal */
void CAssembledRigidModel::getPointOnBodyCurrentCoords(
	const size_t body_index, const mrpt::math::TPoint2D& relative_pt,
	mrpt::math::TPoint2D& out_pt) const
{
	ASSERTDEB_(body_index < m_parent.getBodies().size());

	const CBody& b = m_parent.getBodies()[body_index];

	mrpt::math::TPoint2D q[2];
	this->getPointCurrentCoords(b.points[0], q[0]);
	this->getPointCurrentCoords(b.points[1], q[1]);

	const double L = b.length();
	ASSERTDEB_(L > 0);
	const double Linv = 1.0 / L;

	mrpt::math::TPoint2D u, v;  // unit vectors in X,Y,Z local to the body

	u = (q[1] - q[0]) * Linv;
	v.x = -u.y;
	v.y = u.x;

	out_pt.x = q[0].x + u.x * relative_pt.x + v.x * relative_pt.y;
	out_pt.y = q[0].y + u.y * relative_pt.x + v.y * relative_pt.y;
}

/* Render a ground point */
mrpt::opengl::CSetOfObjects::Ptr
	CAssembledRigidModel::internal_render_ground_point(
		const TMBSPoint& pt, const CBody::TRenderParams& rp) const
{
	mrpt::opengl::CSetOfObjects::Ptr obj =
		mrpt::opengl::CSetOfObjects::Create();
	obj->setLocation(pt.coords.x, pt.coords.y, 0);

	const double support_LXZ = 0.03;
	const double support_LY = 0.05;

	auto gl_box = mrpt::opengl::CBox::Create(
		mrpt::math::TPoint3D(
			-0.5 * support_LXZ, -support_LY, -0.5 * support_LXZ),
		mrpt::math::TPoint3D(0.5 * support_LXZ, 0, 0.5 * support_LXZ), false);
	gl_box->setColor(0, 0, 0.7);
	obj->insert(gl_box);

	return obj;
}

/** Evaluate current energy of the system. */
void CAssembledRigidModel::evaluateEnergy(
	CAssembledRigidModel::TEnergyValues& e) const
{
	timelog.enter("evaluateEnergy");

	e = TEnergyValues();  // Reset to zero

	const std::vector<CBody>& bodies = m_parent.getBodies();
	for (size_t i = 0; i < bodies.size(); i++)
	{
		const CBody& b = bodies[i];

		mrpt::math::TPoint2D dq[2];
		this->getPointCurrentVelocity(b.points[0], dq[0]);
		this->getPointCurrentVelocity(b.points[1], dq[1]);

		const Matrix<double, 2, 1> dq0 = Matrix<double, 2, 1>(dq[0].x, dq[0].y);
		const Matrix<double, 2, 1> dq1 = Matrix<double, 2, 1>(dq[1].x, dq[1].y);

		const Matrix2d& M00 = b.getM00();
		const Matrix2d& M01 = b.getM01();
		const Matrix2d& M11 = b.getM11();

		// Kinetic energy: 0.5 * [q0 q1] * [M00 M10;M10' M11] * [q0 q1]'
		e.E_kin +=
			0.5 * (dq0.transpose() * M00 * dq0 + dq1.transpose() * M11 * dq1)
					  .coeff(0, 0) +
			(dq0.transpose() * M01 * dq1).coeff(0, 0);

		// Potential energy:
		mrpt::math::TPoint2D
			global_cog;  // current COG position, in global coords:
		this->getPointOnBodyCurrentCoords(i, b.cog(), global_cog);

		e.E_pot -= b.mass() * (this->m_gravity[0] * global_cog.x +
							   this->m_gravity[1] * global_cog.y +
							   this->m_gravity[2] * 0 /*global_cog.z*/);
	}

	e.E_total = e.E_kin + e.E_pot;

	timelog.leave("evaluateEnergy");
}
