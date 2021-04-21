/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
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

void ConstraintRelativeAngleAbsolute::update(AssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};
	RelCoordRef angle = actual_rel_coords(arm, 0);

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;

	// Always recalculating L leads to failed numerical Jacobian tests,
	// since it introduces fake dependencies between (x,y) coordinates.
	const double Lsqr_now = Ax * Ax + Ay * Ay;
	if (Lsqr_ == 0 || std::abs(Lsqr_ / Lsqr_now - 1.0) > 0.02)
	{
		// Update cached values
		Lsqr_ = Lsqr_now;
		L_ = std::sqrt(Lsqr_);
	}
	const double L = L_;

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	// This constraints has 2 possible set of equations, with sin() or cos():
	// sin(): pi/4
	// cos(): pi*3/4

	const double theta = angle.x;
	const double w = angle.dotx;
	const double angAcc = angle.ddotx;

	const double sinTh = std::sin(theta), cosTh = std::cos(theta);

	// Use a cached version of "sinTh" to decide which version of the Jacobian
	// to use, since always recalculating it leads to errors in numerical
	// Jacobians:
	if (std::abs(sinTh - sinThCache_) > 0.01)
		// Update cached values
		sinThCache_ = sinTh;
	const bool useCos = std::abs(sinThCache_) > 0.707;

	// Update Phi[i]
	// ----------------------------------
	const double PhiVal = useCos ? Ax - L * cosTh : Ay - L * sinTh;
	arm.Phi_[idx_constr_[0]] = PhiVal;

	// Update dotPhi[i]
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] =
		useCos ? Adotx + L * sinTh * w : Adoty - L * cosTh * w;

	auto& j = jacob.at(0);	// 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (useCos)
	{
		set(j.dPhi_dx[0], -1);
		set(j.dPhi_dx[1], 1);
		set(j.dPhi_drel[0], L * sinTh);

		set(j.dPhi_dy[0], 0);
		set(j.dPhi_dy[1], 0);
	}
	else
	{
		set(j.dPhi_dy[0], -1);
		set(j.dPhi_dy[1], 1);
		set(j.dPhi_drel[0], -L * cosTh);

		set(j.dPhi_dx[0], 0);
		set(j.dPhi_dx[1], 0);
	}

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	set(j.dot_dPhi_dx[0], 0);
	set(j.dot_dPhi_dy[0], 0);
	set(j.dot_dPhi_dx[1], 0);
	set(j.dot_dPhi_dy[1], 0);
	if (useCos)
		set(j.dot_dPhi_drel[0], L * cosTh * w);
	else
		set(j.dot_dPhi_drel[0], L * sinTh * w);

	// Update Phiqq_times_ddq
	// ----------------------------------
	set(j.Phiqq_times_ddq_dx[0], 0);
	set(j.Phiqq_times_ddq_dy[0], 0);
	set(j.Phiqq_times_ddq_dx[1], 0);
	set(j.Phiqq_times_ddq_dy[1], 0);
	if (useCos)
		set(j.Phiqq_times_ddq_drel[0], L * cosTh * angAcc);
	else
		set(j.Phiqq_times_ddq_drel[0], L * sinTh * angAcc);

	// Update dotPhiqq_times_dq
	// ----------------------------------
	set(j.dotPhiqq_times_dq_dx[0], 0);
	set(j.dotPhiqq_times_dq_dy[0], 0);
	set(j.dotPhiqq_times_dq_dx[1], 0);
	set(j.dotPhiqq_times_dq_dy[1], 0);
	if (useCos)
		set(j.dotPhiqq_times_dq_drel[0], -L * w * w * sinTh);
	else
		set(j.dotPhiqq_times_dq_drel[0], L * w * w * cosTh);
}

void ConstraintRelativeAngleAbsolute::print(std::ostream& o) const
{
	o << "ConstraintRelativeAngleAbsolute"
		 "\n";
	o << pointDOFsAsString();
}
