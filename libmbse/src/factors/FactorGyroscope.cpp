/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factors/FactorGyroscope.h>
#include <mbse/AssembledRigidModel.h>

/* #define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif */

using namespace mbse;
using namespace mrpt::math;

gtsam::NonlinearFactor::shared_ptr FactorGyroscope::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorGyroscope::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorGyroscope(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ")\n";
	std::cout << " body: " << body_idx_ << "\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	noiseModel_->print("  noise model: ");
}

bool FactorGyroscope::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol);
}

gtsam::Vector FactorGyroscope::evaluateError(
	const state_t& q_k, const state_t& dq_k, boost::optional<gtsam::Matrix&> H1,
	boost::optional<gtsam::Matrix&> H2) const
{
	const auto n = q_k.size();
	if (dq_k.size() != n)
		throw std::runtime_error("Inconsistent vector lengths!");
	if (n < 1) throw std::runtime_error("Empty state vector!");

	// Set q in the multibody model:
	arm_->q_ = q_k;
	arm_->dotq_ = dq_k;

	const std::vector<Body>& bodies = arm_->mechanism_.bodies();
	ASSERT_BELOW_(body_idx_, bodies.size());

	const Body& body = bodies[body_idx_];

	const size_t pt0_idx = body.points[0];
	const size_t pt1_idx = body.points[1];

	TPoint2D pt0, pt1;
	arm_->getPointCurrentCoords(pt0_idx, pt0);
	arm_->getPointCurrentCoords(pt1_idx, pt1);

	TPoint2D pt0vel, pt1vel;
	arm_->getPointCurrentVelocity(pt0_idx, pt0vel);
	arm_->getPointCurrentVelocity(pt1_idx, pt1vel);

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

	// Evaluate error:
	gtsam::Vector err;
	err.resize(1);

	err[0] = w - reading_;

	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
		Hv.setZero(1, n);

		// point0 & point1:
		const Point2ToDOF pts_dofs[2] = {arm_->points2DOFs_[pt0_idx],
										 arm_->points2DOFs_[pt1_idx]};

		if (size_t i = pts_dofs[0].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt0vel.y - pt1vel.y) -
					   2 * len_inv * len_inv * len_inv * len_inv *
						   (pt0.x - pt1.x) *
						   ((pt0.x - pt1.x) * (pt0vel.y - pt1vel.y) -
							(pt0.y - pt1.y) * (pt0vel.x - pt1vel.x));
		}
		if (size_t i = pts_dofs[0].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = -len_inv * len_inv * (pt0vel.x - pt1vel.x) +
					   2 * len_inv * len_inv * len_inv * len_inv *
						   (pt0.y - pt1.y) *
						   ((pt0.y - pt1.y) * (pt0vel.x - pt1vel.x) -
							(pt0.x - pt1.x) * (pt0vel.y - pt1vel.y));
		}
		if (size_t i = pts_dofs[1].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = -len_inv * len_inv * (pt0vel.y - pt1vel.y) +
					   2 * len_inv * len_inv * len_inv * len_inv *
						   (pt0.x - pt1.x) *
						   ((pt0.x - pt1.x) * (pt0vel.y - pt1vel.y) -
							(pt0.y - pt1.y) * (pt0vel.x - pt1vel.x));
		}
		if (size_t i = pts_dofs[1].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt0vel.x - pt1vel.x) +
					   2 * len_inv * len_inv * len_inv * len_inv *
						   (pt0.y - pt1.y) *
						   (-(pt0.y - pt1.y) * (pt0vel.x - pt1vel.x) +
							(pt0.x - pt1.x) * (pt0vel.y - pt1vel.y));
		}
	}

	// d err / d dq_k
	if (H2)
	{
		auto& Hv = H2.value();
		Hv.setZero(1, n);

		// point0 & point1:
		const Point2ToDOF pts_dofs[2] = {arm_->points2DOFs_[pt0_idx],
										 arm_->points2DOFs_[pt1_idx]};

		if (size_t i = pts_dofs[0].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt1.y - pt0.y);
		}
		if (size_t i = pts_dofs[0].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt0.x - pt1.x);
		}
		if (size_t i = pts_dofs[1].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt0.y - pt1.y);
		}
		if (size_t i = pts_dofs[1].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = len_inv * len_inv * (pt1.x - pt0.x);
		}
	}

	return err;
}
