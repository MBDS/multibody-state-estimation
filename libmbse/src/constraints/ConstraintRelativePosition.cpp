/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/ConstraintRelativePosition.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

/* -------------------------------------------------------------------
				  ConstraintRelativePosition
   -------------------------------------------------------------------*/
void ConstraintRelativePosition::buildSparseStructures(
	AssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	const PointRef p[4] = {
		actual_coords(arm, 0), actual_coords(arm, 1), actual_coords(arm, 2),
		actual_coords(arm, 3)};

	const auto v01 = mrpt::math::TPoint2D(p[1].x - p[0].x, p[1].y - p[0].y);
	const auto v02 = mrpt::math::TPoint2D(p[2].x - p[0].x, p[2].y - p[0].y);
	const auto v03 = mrpt::math::TPoint2D(p[3].x - p[0].x, p[3].y - p[0].y);

	// [ x01 x02 ] [ lambda ] = [ x03 ]
	// [ y01 y02 ] [   mu   ] = [ y03 ]
	// ---- A ----   -- x --  = -- b --

	const Eigen::Matrix2d A =
		(Eigen::Matrix2d() << v01.x, v02.x, v01.y, v02.y).finished();
	const Eigen::Vector2d b = (Eigen::Vector2d() << v03.x, v03.y).finished();

	const Eigen::Vector2d x = A.fullPivLu().solve(b);

	const_cast<ConstraintRelativePosition&>(*this).x_ = x.x();
	const_cast<ConstraintRelativePosition&>(*this).y_ = x.y();
}

void ConstraintRelativePosition::update(AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	const PointRef p[4] = {
		actual_coords(arm, 0), actual_coords(arm, 1), actual_coords(arm, 2),
		actual_coords(arm, 3)};

	const auto v01 = mrpt::math::TPoint2D(p[1].x - p[0].x, p[1].y - p[0].y);
	const auto v02 = mrpt::math::TPoint2D(p[2].x - p[0].x, p[2].y - p[0].y);
	const auto v03 = mrpt::math::TPoint2D(p[3].x - p[0].x, p[3].y - p[0].y);

	const auto dotv01 =
		mrpt::math::TPoint2D(p[1].dotx - p[0].dotx, p[1].doty - p[0].doty);
	const auto dotv02 =
		mrpt::math::TPoint2D(p[2].dotx - p[0].dotx, p[2].doty - p[0].doty);
	const auto dotv03 =
		mrpt::math::TPoint2D(p[3].dotx - p[0].dotx, p[3].doty - p[0].doty);

	// Update Phi[i]
	// ----------------------------------
	arm.Phi_[idx_constr_.at(0)] = v03.x - x_ * v01.x - y_ * v02.x;
	arm.Phi_[idx_constr_.at(1)] = v03.y - x_ * v01.y - y_ * v02.y;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_.at(0)] = dotv03.x - x_ * dotv01.x - y_ * dotv02.x;
	arm.dotPhi_[idx_constr_.at(1)] = dotv03.y - x_ * dotv01.y - y_ * dotv02.y;

	auto& j0 = jacob.at(0);	 // 1st jacob row
	auto& j1 = jacob.at(1);	 // 2nd jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	set(j0.dPhi_dx[0], -1 + x_ + y_);
	set(j0.dPhi_dx[1], -x_);
	set(j0.dPhi_dx[2], -y_);
	set(j0.dPhi_dx[3], 1);

	set(j0.dPhi_dy[0], 0);
	set(j0.dPhi_dy[1], 0);
	set(j0.dPhi_dy[2], 0);
	set(j0.dPhi_dy[3], 0);

	set(j1.dPhi_dy[0], -1 + x_ + y_);
	set(j1.dPhi_dy[1], -x_);
	set(j1.dPhi_dy[2], -y_);
	set(j1.dPhi_dy[3], 1);

	set(j1.dPhi_dx[0], 0);
	set(j1.dPhi_dx[1], 0);
	set(j1.dPhi_dx[2], 0);
	set(j1.dPhi_dx[3], 0);

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	for (int i = 0; i < 4; i++)
	{
		set(j0.dot_dPhi_dx[i], 0);
		set(j0.dot_dPhi_dy[i], 0);
		set(j1.dot_dPhi_dx[i], 0);
		set(j1.dot_dPhi_dy[i], 0);
	}

	// Update Phiqq_times_ddq
	// ----------------------------------
	for (int i = 0; i < 4; i++)
	{
		set(j0.Phiqq_times_ddq_dx[i], 0);
		set(j0.Phiqq_times_ddq_dy[i], 0);
		set(j1.Phiqq_times_ddq_dx[i], 0);
		set(j1.Phiqq_times_ddq_dy[i], 0);
	}

	// Update dotPhiqq_times_dq_dx
	// ----------------------------------
	for (int i = 0; i < 4; i++)
	{
		set(j0.dotPhiqq_times_dq_dx[i], 0);
		set(j0.dotPhiqq_times_dq_dy[i], 0);
		set(j1.dotPhiqq_times_dq_dx[i], 0);
		set(j1.dotPhiqq_times_dq_dy[i], 0);
	}
}

void ConstraintRelativePosition::print(std::ostream& o) const
{
	o << "ConstraintRelativePosition, lambda=" << x_ << ", mu=" << y_ << "\n";
	o << pointDOFsAsString();
}
