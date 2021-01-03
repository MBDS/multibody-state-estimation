/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/CConstraintFixedSlider.h>
#include <mrpt/opengl/CSimpleLine.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

/* -------------------------------------------------------------------
				  CConstraintFixedSlider
   -------------------------------------------------------------------*/
void CConstraintFixedSlider::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	Delta_ = line_pt[1] - line_pt[0];
	ASSERTMSG_(
		Delta_.norm() > 0,
		"Error: the two constraint points must define a line but they "
		"coincide");
	ASSERTMSG_(
		!points_[0]->fixed, "Useless constraint added to a fixed point!");
}

void CConstraintFixedSlider::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates (either fixed or variables in q):
	PointRef p = actual_coords(arm, 0);

	// Update Phi[i] = Ax * (py-y0) - Ay * (px-x0)
	// --------------------------------------------
	const double py_y0 = p.y - line_pt[0].y;
	const double px_x0 = p.x - line_pt[0].x;
	arm.Phi_[idx_constr_[0]] = Delta_.x * py_y0 - Delta_.y * px_x0;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] = Delta_.x * p.doty - Delta_.y * p.dotx;

	auto& j = jacob.at(0);  // 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	set(j.dPhi_dx[0], -Delta_.y);
	set(j.dPhi_dy[0], Delta_.x);

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	if (j.dot_dPhi_dx[0]) *j.dot_dPhi_dx[0] = 0;
	if (j.dot_dPhi_dy[0]) *j.dot_dPhi_dy[0] = 0;

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
