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
#include <mrpt/core/common.h>
#include <mbse/factors/FactorDynamicsIndep.h>
#include <mbse/CAssembledRigidModel.h>
#include <mbse/mbse-utils.h>

#define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorDynamicsIndep::~FactorDynamicsIndep() = default;

gtsam::NonlinearFactor::shared_ptr FactorDynamicsIndep::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorDynamicsIndep::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbde::FactorDynamicsIndep(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ","
			  << keyFormatter(this->key3()) << ")\n";
	noiseModel_->print("  noise model: ");
}

struct NumericJacobParams
{
	CAssembledRigidModel* arm = nullptr;
	CDynamicSimulatorIndepBase* dynamic_solver = nullptr;
	gtsam::Vector z, dz, ddz, q, dq;
};

static void num_err_wrt_z(
	const gtsam::Vector& new_z, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	p.arm->dotq_ = p.dq;
	// Overwrite "z":
	// mbse::insert_subvector(new_z, p.arm->q_, p.arm->get)

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)

	p.dynamic_solver->solve_ddotz(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

static void num_err_wrt_dz(
	const gtsam::Vector& new_dz, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	p.arm->dotq_ = new_dq;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;  // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

bool FactorDynamicsIndep::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   dynamic_solver_->get_model() == e->dynamic_solver_->get_model();
}

gtsam::Vector FactorDynamicsIndep::evaluateError(
	const state_t& z_k, const state_t& dz_k, const state_t& ddz_k,
	boost::optional<gtsam::Matrix&> de_dz,
	boost::optional<gtsam::Matrix&> de_dzp,
	boost::optional<gtsam::Matrix&> de_dzpp) const
{
	MRPT_START

	const auto d = z_k.size();
	ASSERT_EQUAL_(dz_k.size(), z_k.size());
	ASSERT_EQUAL_(ddz_k.size(), z_k.size());
	ASSERT_(z_k.size() > 0);
	ASSERT_(valuesForQk_);

	// Set q & dq in the multibody model:
	CAssembledRigidModel& arm = *dynamic_solver_->get_model_non_const();

	const auto& indepCoordIndices =
		dynamic_solver_->independent_coordinate_indices();
	ASSERT_EQUAL_(indepCoordIndices, z_k.size());

	// Initial guess for "q":
	arm.q_ = valuesForQk_->at<state_t>(key_q_k_);

	// Replace with z and dz:
	arm.set_z(z_k, indepCoordIndices);
	arm.set_dotz(dz_k, indepCoordIndices);
	dynamic_solver_->correct_dependent_q_dq();

	// Predict accelerations:
	Eigen::VectorXd zpp_predicted;
	const double t = 0;  // wallclock time (useless?)
	dynamic_solver_->solve_ddotz(
		t, zpp_predicted, false /* dont auto select indep coords*/);

	// Evaluate error:
	gtsam::Vector err = zpp_predicted - ddz_k;

	// d err / d z_k
	if (de_dz)
	{
		auto& Hv = de_dz.value();
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.dq = dq_k;
		p.ddq = ddq_k;

		const gtsam::Vector x = p.q;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-10);

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
	// d err / d dz_k
	if (de_dzp)
	{
		auto& Hv = de_dzp.value();
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.dq = dq_k;
		p.ddq = ddq_k;

		const gtsam::Vector x = p.dq;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-10);

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
	// d err / d ddz_k
	if (de_dzpp)
	{
		auto& Hv = de_dzpp.value();
		Hv = -Eigen::MatrixXd::Identity(n, n);
	}
	return err;

	MRPT_END
}
