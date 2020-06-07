/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mrpt/core/exceptions.h>
#include <mbse/factors/FactorInverseDynamics.h>
#include <mbse/CAssembledRigidModel.h>

#define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorInverseDynamics::~FactorInverseDynamics() = default;

gtsam::NonlinearFactor::shared_ptr FactorInverseDynamics::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorInverseDynamics::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorInverseDynamics(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ","
			  << keyFormatter(this->key3()) << keyFormatter(this->key4())
			  << ")\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
}

struct NumericJacobParams
{
	CAssembledRigidModel* arm = nullptr;
	CDynamicSimulatorBase* dynamic_solver = nullptr;
	gtsam::Vector q, dq, ddq, Q;
};

static void num_err_wrt_Q(
	const gtsam::Vector& new_Q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->m_q = p.q;
	p.arm->m_dotq = p.dq;
	p.arm->m_Q = new_Q;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

static void num_err_wrt_q(
	const gtsam::Vector& new_q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->m_q = new_q;
	p.arm->m_dotq = p.dq;
	p.arm->m_Q = p.Q;

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
	p.arm->m_Q = p.Q;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

bool FactorInverseDynamics::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   m_dynamic_solver->get_model() == e->m_dynamic_solver->get_model();
}

gtsam::Vector FactorInverseDynamics::evaluateError(
	const state_t& q_k, const state_t& dq_k, const state_t& ddq_k,
	const state_t& Q_k, boost::optional<gtsam::Matrix&> H1,
	boost::optional<gtsam::Matrix&> H2, boost::optional<gtsam::Matrix&> H3,
	boost::optional<gtsam::Matrix&> H4) const
{
	MRPT_START

	const auto n = q_k.size();
	ASSERT_EQUAL_(dq_k.size(), q_k.size());
	ASSERT_EQUAL_(ddq_k.size(), q_k.size());
	ASSERT_EQUAL_(Q_k.size(), q_k.size());
	ASSERT_(q_k.size() > 0);

	// Set q & dq in the multibody model:
	CAssembledRigidModel& arm = *m_dynamic_solver->get_model_non_const();
	arm.m_q = q_k.vector();
	arm.m_dotq = dq_k.vector();
	arm.m_Q = Q_k.vector();

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)

	m_dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	gtsam::Vector err = qpp_predicted - ddq_k.vector();

	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
#if 0
		Hv.setZero(n, n);
#else
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

#endif
	}
	// d err / d dq_k
	if (H2)
	{
		auto& Hv = H2.value();
#if 0
		Hv.setZero(n, n);
#else
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

#endif
	}
	// d err / d ddq_k
	if (H3)
	{
		auto& Hv = H3.value();
		Hv = -Eigen::MatrixXd::Identity(n, n);
	}
	if (H4)
	{
		auto& Hv = H4.value();
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = m_dynamic_solver;
		p.q = q_k.vector();
		p.dq = dq_k.vector();
		p.ddq = ddq_k.vector();
		p.Q = Q_k.vector();

		const gtsam::Vector x = p.Q;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-6);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_Q, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_Q),
			x_incr, p, Hv);
#else
		THROW_EXCEPTION("Implement me!");
#endif
	}
	return err;

	MRPT_END
}
