/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mrpt/core/exceptions.h>
#include <mrpt/core/common.h>
#include <mbse/factors/FactorDynamics.h>
#include <mbse/AssembledRigidModel.h>

#include <gtsam/config.h>
#if defined(GTSAM_USE_TBB)
#error "So far, MBSE is incompatible with GTSAM+TBB!"
#endif
MRPT_TODO(
	"**IMPORTANT** Refactor AssembledRigidModel to separate state and model "
	"data to avoid multithread errors using GTSAM+TBB")

#define USE_NUMERIC_JACOBIAN 1

#if USE_NUMERIC_JACOBIAN
#include <mrpt/math/num_jacobian.h>
#endif

using namespace mbse;

FactorDynamics::~FactorDynamics() = default;

gtsam::NonlinearFactor::shared_ptr FactorDynamics::clone() const
{
	return std::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorDynamics::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorDynamics(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ","
			  << keyFormatter(this->key3()) << ")\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	noiseModel_->print("  noise model: ");
}

struct NumericJacobParams
{
	AssembledRigidModel* arm = nullptr;
	CDynamicSimulatorBase* dynamic_solver = nullptr;
	gtsam::Vector q, dq, ddq;
};

static void num_err_wrt_q(
	const gtsam::Vector& new_q, const NumericJacobParams& p, gtsam::Vector& err)
{
	// Set q & dq in the multibody model:
	p.arm->q_ = new_q;
	p.arm->dotq_ = p.dq;

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

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)
	p.dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	err = qpp_predicted - p.ddq;
}

bool FactorDynamics::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   dynamic_solver_->get_model() == e->dynamic_solver_->get_model();
}

gtsam::Vector FactorDynamics::evaluateError(
	const state_t& q_k, const state_t& dq_k, const state_t& ddq_k,
	gtsam::OptionalMatrixType H1, gtsam::OptionalMatrixType H2,
	gtsam::OptionalMatrixType H3) const
{
	MRPT_START

	const auto n = q_k.size();
	ASSERT_EQUAL_(dq_k.size(), q_k.size());
	ASSERT_EQUAL_(ddq_k.size(), q_k.size());
	ASSERT_(q_k.size() > 0);

	// Set q & dq in the multibody model:
	AssembledRigidModel& arm = *dynamic_solver_->get_model_non_const();
	arm.q_ = q_k;
	arm.dotq_ = dq_k;

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
	const double t = 0;	 // wallclock time (useless?)

	dynamic_solver_->get_model()->realize_operating_point();

	dynamic_solver_->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	gtsam::Vector err = qpp_predicted - ddq_k;

	// d err / d q_k
	if (H1)
	{
		auto& Hv = *H1;
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.dq = dq_k;
		p.ddq = ddq_k;

		const gtsam::Vector x = p.q;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-5);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_q, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_q),
			x_incr, p, Hv);
#else
		Hv.setZero(n, n);
/*
dq = param.dq;
invM = param.invM;
Q = param.Q;
Phiq = ConstLen_q(q);
Phiqq = NumJacob3D(@ConstLen_q, q);
Phiqq_T = NumJacob3D(@ConstLen_T_q, q);
G = Gamma(q,param);
invGamma = inv(G);
c = -InnerTensVect(Phiqq,param.dq)*param.dq;
% ================================== [1]
==================================
A1 = invM*Q;
B1 = Phiq*A1;
C1 = invGamma*B1;
D1 = InnerTensVect(Phiqq_T,C1);
ddq_I = -invM*D1;
dGamma_q = NumJacob3D(@Gamma, q, param);
A2 = InnerTensVect(dGamma_q,C1);
B2 = invGamma*A2;
C2 = Phiq'*B2;
ddq_II = invM*C2;
A3 = invM*Q;
B3 = InnerTensVect(Phiqq,A3);
C3 = invGamma*B3;
D3 = Phiq'*C3;
ddq_III = -invM*D3;
% ================================== [2]
==================================
dotPq = InnerTensVect(Phiqq,dq);A4 = dotPq*dq;
B4 = invGamma*A4;
C4 = InnerTensVect(Phiqq_T,B4);
ddq_IV = -invM*C4;
A5 = InnerTensVect(dGamma_q,B4);
B5 = invGamma*A5;
C5 = Phiq'*B5;
ddq_V = invM*C5;
% Theoretical Jacobian
Jacc_qt = ddq_I+ddq_II+ddq_III+ddq_IV+ddq_V;
*/
#endif
	}
	// d err / d dq_k
	if (H2)
	{
		auto& Hv = *H2;
#if USE_NUMERIC_JACOBIAN
		NumericJacobParams p;
		p.arm = &arm;
		p.dynamic_solver = dynamic_solver_;
		p.q = q_k;
		p.dq = dq_k;
		p.ddq = ddq_k;

		const gtsam::Vector x = p.dq;
		const gtsam::Vector x_incr =
			Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-5);

		mrpt::math::estimateJacobian(
			x,
			std::function<void(
				const gtsam::Vector& new_q, const NumericJacobParams& p,
				gtsam::Vector& err)>(&num_err_wrt_dq),
			x_incr, p, Hv);
#else
		Hv.setZero(n, n);
/*
q = param.q;
invM = param.invM;
Q = param.Q;
Phiq = ConstLen_q(q);
Phiqq = NumJacob3D(@ConstLen_q, q);
Phiqq_T = NumJacob3D(@ConstLen_T_q, q);
G = Gamma(q,param);
invGamma = inv(G);
c = -InnerTensVect(Phiqq,dq)*dq;
dotPq = InnerTensVect(Phiqq,dq);
A = invGamma*dotPq;
B = Phiq'*A;
Jacc_dqt = -2*invM*B;
*/
#endif
	}
	// d err / d ddq_k
	if (H3)
	{
		auto& Hv = *H3;
		Hv = -Eigen::MatrixXd::Identity(n, n);
	}
	return err;

	MRPT_END
}
