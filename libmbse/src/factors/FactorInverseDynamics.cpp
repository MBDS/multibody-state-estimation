/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mrpt/core/exceptions.h>
#include <mbse/factors/FactorInverseDynamics.h>
#include <mbse/AssembledRigidModel.h>

//#define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorInverseDynamics::~FactorInverseDynamics() = default;

gtsam::NonlinearFactor::shared_ptr FactorInverseDynamics::clone() const
{
	return std::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorInverseDynamics::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorInverseDynamics("
			  << keyFormatter(this->key()) << ")\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	// actuatedDegreesOfFreedomInQ
	noiseModel_->print("  noise model: ");
}

#if USE_NUMERIC_JACOBIAN
const double FINITE_DIFF_DELTA = 1e-4;

struct NumericJacobParams
{
	AssembledRigidModel* arm = nullptr;
	CDynamicSimulatorBase* dynamic_solver = nullptr;
	gtsam::Vector q, dq, ddq, Q;
};

static void num_err_wrt_Q(
	const gtsam::Vector& new_Q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	p.arm->dotq_ = p.dq;
	p.arm->Q_ = new_Q;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

static void num_err_wrt_q(
	const gtsam::Vector& new_q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = new_q;
	p.arm->dotq_ = p.dq;
	p.arm->Q_ = p.Q;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

static void num_err_wrt_dq(
	const gtsam::Vector& new_dq, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	p.arm->dotq_ = new_dq;
	p.arm->Q_ = p.Q;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}
#endif

bool FactorInverseDynamics::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   dynamic_solver_->get_model() == e->dynamic_solver_->get_model();
}

gtsam::Vector FactorInverseDynamics::evaluateError(
	const state_t& Q_k, gtsam::OptionalMatrixType d_e_Q) const
{
	MRPT_START

	ASSERT_(valuesFor_q_dq_ != nullptr);

	const auto q_k = valuesFor_q_dq_->at<state_t>(key_q_k_);
	const auto dq_k = valuesFor_q_dq_->at<state_t>(key_dq_k_);
	const auto ddq_k = valuesFor_q_dq_->at<state_t>(key_ddq_k_);

	const auto n = q_k.size();
	ASSERT_EQUAL_(dq_k.size(), q_k.size());
	ASSERT_EQUAL_(ddq_k.size(), q_k.size());
	ASSERT_EQUAL_(Q_k.size(), q_k.size());
	ASSERT_(q_k.size() > 0);

	// Set q & dq in the multibody model:
	// TODO: separate model & state
	AssembledRigidModel& arm = *dynamic_solver_->get_model_non_const();
	arm.q_ = q_k;
	arm.dotq_ = dq_k;
	arm.Q_ = Q_k;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)

	dynamic_solver_->get_model()->realize_operating_point();

	dynamic_solver_->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	gtsam::Vector err = qpp_predicted - ddq_k;

	if (d_e_Q)
	{
		auto& Hv = *d_e_Q;

		const size_t& g = actuatedDegreesOfFreedomInQ_.size();

		Eigen::MatrixXd R(n, g);
		for (size_t i = 0; i < g; i++)
		{
			AssembledRigidModel::ComputeDependentParams dp;
			AssembledRigidModel::ComputeDependentResults dr;

			arm.dotq_.setZero();
			arm.dotq_[actuatedDegreesOfFreedomInQ_[i]] = 1.0;
			arm.computeDependentPosVelAcc(
				actuatedDegreesOfFreedomInQ_, false, true, dp, dr);
			R.col(i) = arm.dotq_;
		}
		arm.dotq_ = dq_k;

		const Eigen::MatrixXd M = arm.buildMassMatrix_dense();
		const auto Mr = R.transpose() * M * R;

		const Eigen::MatrixXd J_z = Mr.inverse() * R.transpose();

		Hv.setZero(n, n);
		for (size_t i = 0; i < g; i++)
		{
			Hv.row(actuatedDegreesOfFreedomInQ_[i]) = J_z.row(i);
		}

		// std::cout << "\nde/d{Q} THEORETHICAL:\n" << Hv << "\n\n";

#if USE_NUMERIC_JACOBIAN

		const double cacheTol_q = 0.01;
		const double cacheTol_dq = 0.1;
		const double cacheTol_ddq = 0.5;

		if (cached_q_.size() == q_k.size() &&
			(cached_q_ - q_k).array().abs().maxCoeff() < cacheTol_q &&
			(cached_dq_ - dq_k).array().abs().maxCoeff() < cacheTol_dq &&
			(cached_ddq_ - ddq_k).array().abs().maxCoeff() < cacheTol_ddq)
		{
			Hv = cached_d_e_Q_;
		}
		else
		{
			NumericJacobParams p;
			p.arm = &arm;
			p.dynamic_solver = dynamic_solver_;
			p.q = q_k;
			p.dq = dq_k;
			p.ddq = ddq_k;
			p.Q = Q_k;

			const gtsam::Vector x = p.Q;
			const gtsam::Vector x_incr = Eigen::VectorXd::Constant(
				x.rows(), x.cols(), FINITE_DIFF_DELTA);

			mrpt::math::estimateJacobian(
				x,
				std::function<void(
					const gtsam::Vector& new_Q, const NumericJacobParams& p,
					gtsam::Vector& err)>(&num_err_wrt_Q),
				x_incr, p, Hv);

			cached_d_e_Q_ = Hv;
			cached_q_ = q_k;
			cached_dq_ = dq_k;
			cached_ddq_ = ddq_k;
		}

		std::cout << "\nde/d{Q} NUMERIC:\n" << Hv << "\n\n";
#endif
	}
	return err;

	MRPT_END
}
