/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CModelDefinition.h>
#include <mbse/CAssembledRigidModel.h>
#include <mrpt/core/exceptions.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

void CAssembledRigidModel::builGeneralizedForces(Eigen::VectorXd& Q) const
{
	const size_t nDOFs = m_q.size();
	Q.resize(nDOFs);
	builGeneralizedForces(&Q[0]);
}

/* -------------------------------------------------------------------
				  builGeneralizedForces
-------------------------------------------------------------------*/
void CAssembledRigidModel::builGeneralizedForces(double* q) const
{
	timelog.enter("builGeneralizedForces");

	const size_t nDOFs = m_q.size();
	Eigen::Map<Eigen::VectorXd> Q(q, nDOFs);
	Q.setZero();

	const std::vector<CBody>& parent_bodies = m_parent.getBodies();

	// Gravity force:
	// --------------------------------
	// For each body:
	for (size_t i = 0; i < parent_bodies.size(); i++)
	{
		const CBody& body = parent_bodies[i];

		// Compute the body local coordinates (a,b) of the point at
		// where the force is applied:
		TPoint2D force_local_point;

		MRPT_TODO("Get current coordinates, from m_q!")
		const TMBSPoint& p0_info = m_parent.getPointInfo(body.points[0]);
		const TMBSPoint& p1_info = m_parent.getPointInfo(body.points[1]);

		const TPoint2DOF& p0_dofs = m_points2DOFs[body.points[0]];
		const TPoint2DOF& p1_dofs = m_points2DOFs[body.points[1]];

		if (0)
		{
			// Get the current global coordinates of the body points:

			const TPoint2D p0(
				(p0_dofs.dof_x != INVALID_DOF) ? m_q[p0_dofs.dof_x]
											   : p0_info.coords.x,
				(p0_dofs.dof_y != INVALID_DOF) ? m_q[p0_dofs.dof_y]
											   : p0_info.coords.y);
			const TPoint2D p1(
				(p1_dofs.dof_x != INVALID_DOF) ? m_q[p1_dofs.dof_x]
											   : p1_info.coords.x,
				(p1_dofs.dof_y != INVALID_DOF) ? m_q[p1_dofs.dof_y]
											   : p1_info.coords.y);

			// const TPoint2D p1...
		}
		else
		{
			// Gravity force is always applied at the cog, whose coordinates are
			// ALREADY stored as LOCAL COORDINATES:
			force_local_point = body.cog();
		}

		// Cp matrix. Eq. (62), pag. 105 from J. Cuadrado's manual.
		const double a = force_local_point.x;
		const double b = force_local_point.y;
		const double L = body.length();

		const double Cp_vals[2 * 4] = {L - a, b, a, -b, -b, L - a, b, a};
		Eigen::Matrix<double, 2, 4, Eigen::RowMajor> Cp(Cp_vals);
		Cp *= 1.0 / L;

		// Force vector:
		Eigen::Vector2d F;
		// F = body.mass * m_gravity;
		F[0] = body.mass() * m_gravity[0];	// x
		F[1] = body.mass() * m_gravity[1];	// y

		// Q = Cp^t * F
		const Eigen::Vector4d Qi = Cp.transpose() * F;

		// Assemble:
		const bool p0_fixed = p0_info.fixed;
		const bool p1_fixed = p1_info.fixed;

		const size_t idx_x0 =
			p0_dofs.dof_x;	// Will be INVALID_DOF if it's a fixed point
		const size_t idx_x1 = p1_dofs.dof_x;

		if (!p0_fixed) Q.segment<2>(idx_x0) += Qi.head<2>();
		if (!p1_fixed) Q.segment<2>(idx_x1) += Qi.tail<2>();
	}

	// External forces:
	// --------------------------------
	ASSERT_EQUAL_(Q.rows(), m_Q.rows());
	ASSERT_EQUAL_(Q.cols(), m_Q.cols());
	Q += m_Q;

	timelog.leave("builGeneralizedForces");
}
