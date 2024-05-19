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
#include <mrpt/core/common.h>
#include <mbse/factors/FactorDynamicsIndep.h>
#include <mbse/AssembledRigidModel.h>
#include <mbse/mbse-utils.h>

#define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorDynamicsIndep::~FactorDynamicsIndep() = default;

gtsam::NonlinearFactor::shared_ptr FactorDynamicsIndep::clone() const
{
return gtsam::NonlinearFactor::shared_ptr(new This(*this));
}

void FactorDynamicsIndep::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorDynamicsIndep(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ","
			  << keyFormatter(this->key3()) << ")\n";
	noiseModel_->print("  noise model: ");
}

struct NumericJacobParams
{
	AssembledRigidModel* arm = nullptr;
	CDynamicSimulatorIndepBase* dynamic_solver = nullptr;
	gtsam::Vector z, dz, ddz, q;
};

static void num_err_wrt_z(
	const gtsam::Vector& new_z, const NumericJacobParams& p, gtsam::Vector& err)
{
	const auto& zIndices = p.dynamic_solver->independent_coordinate_indices();

	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	// Overwrite "z" & "dz":
	mbse::overwrite_subset(p.arm->q_, new_z, zIndices);
	mbse::overwrite_subset(p.arm->dotq_, p.dz, zIndices);

	// Ensure q and dq are updated after the change in "z":
	AssembledRigidModel::ComputeDependentParams cdp;
	AssembledRigidModel::ComputeDependentResults cdr;
	cdp.nItersMax = 3;
	p.arm->computeDependentPosVelAcc(
		zIndices, true /*update_q*/, true /* update_dq*/, cdp, cdr);

	ASSERT_LT_(cdr.pos_final_phi, 1e-3);

	// Predict accelerations:
	const double t = 0;	 // wallclock time (useless?)
	Eigen::VectorXd zpp_predicted;

	p.dynamic_solver->solve_ddotz(t, zpp_predicted);

	// Evaluate error:
	err = zpp_predicted - p.ddz;
}

static void num_err_wrt_dz(
	const gtsam::Vector& new_dz, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	const auto& zIndices = p.dynamic_solver->independent_coordinate_indices();

	// Set q & dq in the multibody model:
	p.arm->q_ = p.q;
	// Overwrite "z" & "dz":
	mbse::overwrite_subset(p.arm->q_, p.z, zIndices);
	mbse::overwrite_subset(p.arm->dotq_, new_dz, zIndices);

	// Ensure q and dq are updated after the change in "z":
	AssembledRigidModel::ComputeDependentParams cdp;
	AssembledRigidModel::ComputeDependentResults cdr;
	cdp.nItersMax = 3;
	p.arm->computeDependentPosVelAcc(
		zIndices, true /*update_q*/, true /* update_dq*/, cdp, cdr);

	ASSERT_LT_(cdr.pos_final_phi, 1e-3);

	// Predict accelerations:
	const double t = 0;	 // wallclock time (useless?)
	Eigen::VectorXd zpp_predicted;

	p.dynamic_solver->solve_ddotz(t, zpp_predicted);

	// Evaluate error:
	err = zpp_predicted - p.ddz;
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
	boost::optional<gtsam::Matrix&>  de_dz, boost::optional<gtsam::Matrix&>  de_dzp,
	boost::optional<gtsam::Matrix&>  de_dzpp) const
{
	MRPT_START

	const auto d = z_k.size();
	ASSERT_EQUAL_(dz_k.size(), z_k.size());
	ASSERT_EQUAL_(ddz_k.size(), z_k.size());
	ASSERT_(z_k.size() > 0);
	ASSERT_(valuesForQk_);

	// Set q & dq in the multibody model:
	AssembledRigidModel& arm = *dynamic_solver_->get_model_non_const();

	const auto& indepCoordIndices =
		dynamic_solver_->independent_coordinate_indices();
	ASSERT_EQUAL_(
		static_cast<size_t>(indepCoordIndices.size()),
		static_cast<size_t>(z_k.size()));

	// Initial guess for "q":
	const Eigen::VectorXd q_k = valuesForQk_->at<state_t>(key_q_k_);
	arm.q_ = q_k;

	// Replace with z and dz:
	mbse::overwrite_subset(arm.q_, z_k, indepCoordIndices);
	mbse::overwrite_subset(arm.dotq_, dz_k, indepCoordIndices);

	const double fdErr = arm.finiteDisplacement(
		indepCoordIndices, 1e-9, 10 /*max iters*/,
		true /* also solve dot{q} */);
	ASSERT_LT_(fdErr, 1e-3);

	// Predict accelerations:
	Eigen::VectorXd zpp_predicted;
	const double t = 0;	 // wallclock time (useless?)

	dynamic_solver_->get_model()->realize_operating_point();
	dynamic_solver_->solve_ddotz(t, zpp_predicted);

	// Evaluate error:
	gtsam::Vector err = zpp_predicted - ddz_k;

	// d err / d z_k
	if (de_dz)
	{
		auto& Hv = *de_dz;
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.z = z_k;
		p.dz = dz_k;
		p.ddz = ddz_k;

		const gtsam::Vector x = p.z;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-5);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_q, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_z),
			x_incr, p, Hv);
#else
		Hv.setZero(n, n);
#endif
	}
	// d err / d dz_k
	if (de_dzp)
	{
		auto& Hv = *de_dzp;
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.z = z_k;
		p.dz = dz_k;
		p.ddz = ddz_k;

		const gtsam::Vector x = p.dz;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-5);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_q, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_dz),
			x_incr, p, Hv);
#else
		Hv.setZero(n, n);
#endif
	}
	// d err / d ddz_k
	if (de_dzpp)
	{
		auto& Hv = *de_dzpp;
		Hv = -Eigen::MatrixXd::Identity(d, d);
	}
	return err;

	MRPT_END
}
