/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factors/FactorConstraintsVel.h>
#include <mbse/CAssembledRigidModel.h>
#include <mrpt/core/exceptions.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#define USE_NUMERIC_JACOBIAN 0

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorConstraintsVel::~FactorConstraintsVel() = default;

gtsam::NonlinearFactor::shared_ptr FactorConstraintsVel::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorConstraintsVel::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorConstraintsVel("
			  << keyFormatter(this->key1()) << "," << keyFormatter(this->key2())
			  << ")\n";
	noiseModel_->print("  noise model: ");
}

bool FactorConstraintsVel::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol);
}

#if USE_NUMERIC_JACOBIAN
struct NumericJacobParams
{
	CAssembledRigidModel* arm = nullptr;
	gtsam::Vector q, dq;
};

static void num_err_wrt_q(
	const gtsam::Vector& new_q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = new_q;
	p.arm->dotq_ = p.dq;

	p.arm->update_numeric_Phi_and_Jacobians();

	// Evaluate error:
	const Eigen::MatrixXd Phi_q = p.arm->Phi_q_.asDense();
	err = Phi_q * p.arm->dotq_;
}

static void num_err_wrt_dq(
	const gtsam::Vector& new_dq, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	p.arm->dotq_ = new_dq;

	p.arm->update_numeric_Phi_and_Jacobians();

	// Evaluate error:
	const Eigen::MatrixXd Phi_q = p.arm->Phi_q_.asDense();
	err = Phi_q * p.arm->dotq_;
}
#endif

gtsam::Vector FactorConstraintsVel::evaluateError(
	const state_t& q_k, const state_t& dotq_k,
	boost::optional<gtsam::Matrix&> H1,
	boost::optional<gtsam::Matrix&> H2) const
{
	MRPT_START

	// const auto n = q_k.size();
	// const auto m = arm_->Phi_.rows();

	ASSERT_EQUAL_(dotq_k.size(), q_k.size());
	ASSERT_(q_k.size() > 0);

	// Set q in the multibody model:
	arm_->q_ = q_k;
	arm_->dotq_ = dotq_k;

	// Update Jacobian and Hessian tensor:
	arm_->update_numeric_Phi_and_Jacobians();

	// Evaluate error:
	const Eigen::MatrixXd Phi_q = arm_->Phi_q_.asDense();
	const Eigen::MatrixXd dPhiqdq_dq = arm_->dotPhi_q_.asDense();

	gtsam::Vector err = Phi_q * dotq_k;

	// Get the Jacobians required for optimization:
	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = arm_.get();
		p.q = q_k;
		p.dq = dotq_k;

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
		Hv = dPhiqdq_dq;

#endif
	}

	if (H2)
	{
		auto& Hv = H2.value();
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = arm_.get();
		p.q = q_k;
		p.dq = dotq_k;

		const gtsam::Vector x = p.dq;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-6);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_dq, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_dq),
			x_incr, p, Hv);
#else
		Hv = Phi_q;
#endif
	}

	return err;

	MRPT_END
}
