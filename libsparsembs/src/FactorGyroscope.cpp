#include <sparsembs/FactorGyroscope.h>
#include <sparsembs/CAssembledRigidModel.h>

/* #define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif */

using namespace sparsembs;

gtsam::NonlinearFactor::shared_ptr FactorGyroscope::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorGyroscope::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorGyroscope(" << keyFormatter(this->key1()) << ","
			  << keyFormatter(this->key2()) << ")\n";
	std::cout << " body: " << m_body_idx << "\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
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
	m_arm->m_q = q_k;
	m_arm->m_dotq = dq_k;

	const std::vector<CBody>& bodies = m_arm->m_parent.getBodies();
	ASSERT_BELOW_(m_body_idx, bodies.size());

	const CBody& body = bodies[m_body_idx];

	const size_t pt0_idx = body.points[0];
	const size_t pt1_idx = body.points[1];

	TPoint2D pt0, pt1;
	m_arm->getPointCurrentCoords(pt0_idx, pt0);
	m_arm->getPointCurrentCoords(pt1_idx, pt1);

	TPoint2D pt0vel, pt1vel;
	m_arm->getPointCurrentVelocity(pt0_idx, pt0vel);
	m_arm->getPointCurrentVelocity(pt1_idx, pt1vel);

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

	err[0] = w - m_reading;

	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
		Hv.setZero(1, n);

		// point0 & point1:
		const TMBSPoint pts_info[2] = {m_arm->m_parent.getPointInfo(pt0_idx),
									   m_arm->m_parent.getPointInfo(pt1_idx)};
		const TPoint2DOF pts_dofs[2] = {m_arm->m_points2DOFs[pt0_idx],
										m_arm->m_points2DOFs[pt1_idx]};

		if (size_t i = pts_dofs[0].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt0vel.y - pt1vel.y;  // dq_k[1] - dq_k[3];
		}
		if (size_t i = pts_dofs[0].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt1vel.x - pt0vel.x;  // dq_k[2] - dq_k[0];
		}
		if (size_t i = pts_dofs[1].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt1vel.y - pt0vel.y;  // dq_k[3] - dq_k[1];
		}
		if (size_t i = pts_dofs[1].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt1vel.x - pt0vel.x;  // dq_k[2] - dq_k[0];
		}
	}

	// d err / d dq_k
	if (H2)
	{
		auto& Hv = H2.value();
		Hv.setZero(1, n);

		// point0 & point1:
		const TMBSPoint pts_info[2] = {m_arm->m_parent.getPointInfo(pt0_idx),
									   m_arm->m_parent.getPointInfo(pt1_idx)};
		const TPoint2DOF pts_dofs[2] = {m_arm->m_points2DOFs[pt0_idx],
										m_arm->m_points2DOFs[pt1_idx]};

		if (size_t i = pts_dofs[0].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt1.y - pt0.y;  // q_k[3] - q_k[1];
		}
		if (size_t i = pts_dofs[0].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt0.x - pt1.x;  // q_k[0] - q_k[2];
		}
		if (size_t i = pts_dofs[1].dof_x; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt0.y - pt1.y;  // q_k[1] - q_k[3];
		}
		if (size_t i = pts_dofs[1].dof_y; i != INVALID_DOF)
		{
			// x0 is NOT a fixed point, it DO belong to q:
			// fill Jacobian for column "i"
			Hv(0, i) = pt1.x - pt0.x;  // q_k[2] - q_k[0];
		}
	}

	return err;
}
