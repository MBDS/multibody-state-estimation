/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/AssembledRigidModel.h>
#include <mbse/virtual-sensors.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------
// Virtual sensor: Gyroscope
// ---------------------------------------
double CVirtualSensor_Gyro::simulate_reading(
	const AssembledRigidModel& arm) const
{
	// Estimate the current angular velocity of the i'th body:
	const std::vector<Body>& bodies = arm.mechanism_.bodies();
	ASSERTDEB_(body_idx_ < bodies.size());

	const Body& body = bodies[body_idx_];

	const size_t pt0_idx = body.points[0];
	const size_t pt1_idx = body.points[1];

	TPoint2D pt0, pt1;
	arm.getPointCurrentCoords(pt0_idx, pt0);
	arm.getPointCurrentCoords(pt1_idx, pt1);

	TPoint2D pt0vel, pt1vel;
	arm.getPointCurrentVelocity(pt0_idx, pt0vel);
	arm.getPointCurrentVelocity(pt1_idx, pt1vel);

	// u: unit director vector from pt0->pt1
	TPoint2D u = pt1 - pt0;
	const double len = u.norm();
	const double len_inv = 1.0 / len;
	u *= len_inv;

	// +90deg orthogonal direction:
	const TPoint2D v(-u.y, u.x);

	// The relative velocity of pt1 wrt pt0:
	const TPoint2D rel_vel = pt1vel - pt0vel;

	// Its one-dimensional value is the projection of that vector on "v":
	const double rel_vel_value = rel_vel.x * v.x + rel_vel.y * v.y;

	// w = rel_vel / distance
	const double w = rel_vel_value * len_inv;

	return w;
}
