/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/constraints/CConstraintMobileSlider.h>

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
	// Get references to the point coordinates and velocities
	// (either fixed or variables in q):
	PointRef p = actual_coords(arm, 0);
	PointRef pr[2] = {actual_coords(arm, 1), actual_coords(arm, 2)};

	// Update Phi[i]
	// ----------------------------------
	arm.Phi_[idx_constr_[0]] = (pr[1].x - pr[0].x) * (p.y - pr[0].y) -
							   (pr[1].y - pr[0].y) * (p.x - pr[0].x);

	// Update dotPhi[i] (partial-Phi[i]_partial-t)
	// ----------------------------------
	arm.dotPhi_[idx_constr_[0]] = (pr[1].dotx - pr[0].dotx) * (p.y - pr[0].y) +
								  (pr[1].x - pr[0].x) * (p.doty - pr[0].doty) -
								  (pr[1].doty - pr[0].doty) * (p.x - pr[0].x) -
								  (pr[1].y - pr[0].y) * (p.dotx - pr[0].dotx);

	auto& j = jacob.at(0);	// 1st (and unique) jacob row

	// Update Jacobian dPhi_dq(i,:)
	// ----------------------------------
	if (j.dPhi_dx[0]) *j.dPhi_dx[0] = pr[0].y - pr[1].y;
	if (j.dPhi_dy[0]) *j.dPhi_dy[0] = pr[1].x - pr[0].x;

	if (j.dPhi_dx[1]) *j.dPhi_dx[1] = -p.y + pr[1].y;
	if (j.dPhi_dy[1]) *j.dPhi_dy[1] = p.x - pr[1].x;

	if (j.dPhi_dx[2]) *j.dPhi_dx[2] = p.y - pr[0].y;
	if (j.dPhi_dy[2]) *j.dPhi_dy[2] = -p.x + pr[0].x;

	// Update Jacobian \dot{dPhi_dq}(i,:)
	// ----------------------------------
	set(j.dot_dPhi_dx[0], pr[0].doty - pr[1].doty);
	set(j.dot_dPhi_dy[0], pr[1].dotx - pr[0].dotx);

	set(j.dot_dPhi_dx[1], -p.doty + pr[1].doty);
	set(j.dot_dPhi_dy[1], p.dotx - pr[1].dotx);

	set(j.dot_dPhi_dx[2], p.doty - pr[0].doty);
	set(j.dot_dPhi_dy[2], -p.dotx + pr[0].dotx);
}
