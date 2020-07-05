/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/CConstraintRelativeAngleAbsolute.h>
#include <mrpt/opengl/CSimpleLine.h>
#include <mrpt/math/wrap2pi.h>

using namespace mbse;
using namespace Eigen;
using mrpt::square;

void CConstraintRelativeAngleAbsolute::buildSparseStructures(
	CAssembledRigidModel& arm) const
{
	commonbuildSparseStructures(arm);

	ASSERTMSG_(
		!(points_[0]->fixed && points_[1]->fixed),
		"Useless relative coordinate added between two fixed points!");
}

void CConstraintRelativeAngleAbsolute::update(CAssembledRigidModel& arm) const
{
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p[2] = {actual_coords(arm, 0), actual_coords(arm, 1)};
	RelCoordRef angle = actual_rel_coords(arm, 0);

	const double Ax = p[1].x - p[0].x;
	const double Ay = p[1].y - p[0].y;
	const double L = std::sqrt(Ax * Ax + Ay * Ay);

	const double Adotx = p[1].dotx - p[0].dotx;
	const double Adoty = p[1].doty - p[0].doty;

	// This constraints has 2 possible set of equations, with sin() or cos():
	// sin(): pi/4
	// cos(): pi*3/4

	const double theta = angle.x;
	const double w = angle.dotx;

	const double sinTh = std::sin(theta), cosTh = std::cos(theta);

	const bool useCos = std::abs(sinTh) > 0.707;

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

	// Update Jacobian \{dPhiq*dq}_{\dot{q}}(i,:)
	// --------------------------------------
	if (useCos)
	{
		set(j.dPhiqdq_dx[0], -1);
		set(j.dPhiqdq_dx[1], 1);
		set(j.dPhiqdq_dy[0], 0);
		set(j.dPhiqdq_dy[1], 0);
		set(j.dPhiqdq_drel[0], L * cosTh);
	}
	else
	{
		set(j.dPhiqdq_dx[0], 0);
		set(j.dPhiqdq_dx[1], 0);
		set(j.dPhiqdq_dy[0], -1);
		set(j.dPhiqdq_dy[1], 1);
		set(j.dPhiqdq_drel[0], L * sinTh);
	}
}
