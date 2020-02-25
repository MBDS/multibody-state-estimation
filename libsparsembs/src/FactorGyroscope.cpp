#include <sparsembs/FactorGyroscope.h>
#include <sparsembs/CAssembledRigidModel.h>

/* #define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif */

using namespace sparsembs;

FactorGyroscope::~FactorGyroscope() = default;

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
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
}

/*  struct NumericJacobParams
{
	CAssembledRigidModel* arm = nullptr;
	CDynamicSimulatorBase* dynamic_solver = nullptr;
	gtsam::Vector q, dq, ddq;
};

static void num_err_wrt_q(
	const gtsam::Vector& new_q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->m_q = new_q;
	p.arm->m_dotq = p.dq;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

static void num_err_wrt_dq(
	const gtsam::Vector& new_dq, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->m_q = p.q;
	p.arm->m_dotq = new_dq;

	// Predict angular velocity:
	Eigen::VectorXd omega_predicted;
	const double t = 0;  // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
} */

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

	// Set q & dq in the multibody model:
	CAssembledRigidModel& arm =
		*m_gyro->simulate_reading(const CAssembledRigidModel& dq_k);
	/* arm.m_q = q_k.vector();
	arm.m_dotq = dq_k.vector(); */
	/*
		// Predict angular velocity:
		Eigen::VectorXd qp_predicted;
		const double t = 0;  // wallclock time (useless?)

		m_dynamic_solver->solve_ddotq(t, qpp_predicted);

		// Evaluate error:
		gtsam::Vector err = qpp_predicted - ddq_k.vector();

		// d err / d q_k
		if (H1)
		{
			auto& Hv = H1.value();
	#if USE_NUMERIC_JACOBIAN
			NumericJacobParams p;
			p.arm = &arm;
			p.dynamic_solver = m_dynamic_solver;
			p.q = q_k.vector();
			p.dq = dq_k.vector();
			p.ddq = ddq_k.vector();

			const gtsam::Vector x = p.q;
			const gtsam::Vector x_incr =
				Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-6);

			mrpt::math::estimateJacobian(
				x,
				std::function<void(
					const gtsam::Vector& new_q, const NumericJacobParams& p,
					gtsam::Vector& err)>(&num_err_wrt_q),
				x_incr, p, Hv);
	#else
			Hv.setZero(n, n);
	#endif
		}
		// d err / d dq_k
		if (H2)
		{
			auto& Hv = H2.value();
	/* #if USE_NUMERIC_JACOBIAN
			NumericJacobParams p;
			p.arm = &arm;
			p.dynamic_solver = m_dynamic_solver;
			p.q = q_k.vector();
			p.dq = dq_k.vector();
			p.ddq = ddq_k.vector();

			const gtsam::Vector x = p.dq;
			const gtsam::Vector x_incr =
				Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-6);

			mrpt::math::estimateJacobian(
				x,
				std::function<void(
					const gtsam::Vector& new_q, const NumericJacobParams& p,
					gtsam::Vector& err)>(&num_err_wrt_dq),
				x_incr, p, Hv);
	#else
			Hv.setZero(n, n);
	#endif
		}
		// d err / d ddq_k
		if (H3)
		{
			auto& Hv = H3.value();
			Hv = -Eigen::MatrixXd::Identity(n, n);
		}
		return err;
	} */
