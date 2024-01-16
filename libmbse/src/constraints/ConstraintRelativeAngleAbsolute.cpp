/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/ConstraintRelativeAngleAbsolute.h>
#include <mrpt/opengl/CSimpleLine.h>
#include <mrpt/math/wrap2pi.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

void ConstraintRelativeAngleAbsolute::buildSparseStructures(
	AssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless relative coordinate added between two fixed points!");
}

void ConstraintRelativeAngleAbsolute::realizeOperatingPoint(
	const AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};
	RelCoordRef angle = actual_rel_coords(arm, 0);

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	// Always recalculating L leads to failed numerical Jacobian tests,
	// since it introduces fake dependencies between (x,y) coordinates.
	if (L_ == 0)
	{
		const auto Lsqr = Ax * Ax + Ay * Ay;
		L_ = std::sqrt(Lsqr);
		// std::cout << "UPDATING L: " << L_ << std::endl;
	}

	const double theta = angle.x;
	const double sinTh = std::sin(theta);

	useCos_ = std::abs(sinTh) > 0.707;
}

void ConstraintRelativeAngleAbsolute::update(AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};
	RelCoordRef angle = actual_rel_coords(arm, 0);

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	// This constraints has 2 possible set of equations, with sin() or cos():
	// sin(): pi/4
	// cos(): pi*3/4

	const double theta = angle.x;
	const double w = angle.dotx;
	const double angAcc = angle.ddotx;
	const double sinTh = std::sin(theta), cosTh = std::cos(theta);

	// Update Phi[i]
	// ----------------------------------
	const double PhiVal = useCos_ ? Ax - L_ * cosTh : Ay - L_ * sinTh;
	arm.Phi_[idx_constr_[0]] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] =
		useCos_ ? Adotx + L_ * sinTh * w : Adoty - L_ * cosTh * w;

	auto& j = jacob.at(0);	// 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (useCos_)
	{
		set(j.dPhi_dx[0], -1);
		set(j.dPhi_dx[1], 1);
		set(j.dPhi_drel[0], L_ * sinTh);

		set(j.dPhi_dy[0], 0);
		set(j.dPhi_dy[1], 0);
	}
	else
	{
		set(j.dPhi_dy[0], -1);
		set(j.dPhi_dy[1], 1);
		set(j.dPhi_drel[0], -L_ * cosTh);

		set(j.dPhi_dx[0], 0);
		set(j.dPhi_dx[1], 0);
	}

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	set(j.dot_dPhi_dx[0], 0);
	set(j.dot_dPhi_dy[0], 0);
	set(j.dot_dPhi_dx[1], 0);
	set(j.dot_dPhi_dy[1], 0);
	if (useCos_)
		set(j.dot_dPhi_drel[0], L_ * cosTh * w);
	else
		set(j.dot_dPhi_drel[0], L_ * sinTh * w);

	// Update Phiqq_times_ddq
	// ----------------------------------
	set(j.Phiqq_times_ddq_dx[0], 0);
	set(j.Phiqq_times_ddq_dy[0], 0);
	set(j.Phiqq_times_ddq_dx[1], 0);
	set(j.Phiqq_times_ddq_dy[1], 0);
	if (useCos_)
		set(j.Phiqq_times_ddq_drel[0], L_ * cosTh * angAcc);
	else
		set(j.Phiqq_times_ddq_drel[0], L_ * sinTh * angAcc);

	// Update dotPhiqq_times_dq
	// ----------------------------------
	set(j.dotPhiqq_times_dq_dx[0], 0);
	set(j.dotPhiqq_times_dq_dy[0], 0);
	set(j.dotPhiqq_times_dq_dx[1], 0);
	set(j.dotPhiqq_times_dq_dy[1], 0);
	if (useCos_)
		set(j.dotPhiqq_times_dq_drel[0], -L_ * w * w * sinTh);
	else
		set(j.dotPhiqq_times_dq_drel[0], L_ * w * w * cosTh);
}

void ConstraintRelativeAngleAbsolute::print(std::ostream& o) const
{
	o << "ConstraintRelativeAngleAbsolute"
		 "\n";
	o << pointDOFsAsString();
}
