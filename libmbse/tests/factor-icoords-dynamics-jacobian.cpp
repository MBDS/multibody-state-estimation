/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <gtest/gtest.h>

#include <mbse/model-examples.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/CAssembledRigidModel.h>
#include <gtsam/inference/Symbol.h>
#include <mbse/factors/FactorDynamicsIndep.h>
#include <mrpt/math/num_jacobian.h>
#include <mrpt/system/CTimeLogger.h>

using namespace std;
using namespace mbse;

struct NumericJacobParams
{
	gtsam::Vector q0, z, dz, ddz;

	/** Especifies which variable are we taking numerical derivatives with
	 * respect to: 0: q, 1: \dot{q}, 2: \ddot{q}
	 */
	int diff_variable = 0;
	gtsam::NoiseModelFactor3<state_t, state_t, state_t>* factor = nullptr;
};

static void num_err_wrt_state(
	const gtsam::Vector& new_state, const NumericJacobParams& p,
	gtsam::Vector& err)
{
	auto z = state_t(p.diff_variable == 0 ? new_state : p.z);
	auto dz = state_t(p.diff_variable == 1 ? new_state : p.dz);
	auto ddz = state_t(p.diff_variable == 2 ? new_state : p.ddz);

	// Evaluate error:
	err = p.factor->evaluateError(z, dz, ddz);
}

auto alternative_symbol(const char c)
{
	return [c](std::uint64_t j) { return gtsam::Symbol(c, j); };
}

TEST(Jacobians, FactorDynamicsIndepCoords)
{
	try
	{
		const auto DDZ = alternative_symbol('a');
		const auto DZ = alternative_symbol('v');
		const auto Z = alternative_symbol('q');
		const auto Q = gtsam::symbol_shorthand::Q;

		using namespace mbse;

		// Create the multibody object:
		CModelDefinition model;
		std::vector<RelativeDOF> rDOFs;

		mbse::buildFourBarsMBS(model);
		// Add an extra relative coordinate:
		rDOFs.emplace_back(mbse::RelativeAngleAbsoluteDOF(0, 1));

		std::shared_ptr<CAssembledRigidModel> aMBS =
			model.assembleRigidMBS(rDOFs);
		aMBS->setGravityVector(0, -9.81, 0);

		aMBS->printCoordinates();

		CDynamicSimulator_Indep_dense dynSimul(aMBS);

		// Enforce the use of the theta angle as independent coordinate:
		dynSimul.independent_coordinate_indices({4});

		// Must be called before solve_ddotq():
		dynSimul.prepare();

		// Add factors:
		// Create factor noises:
		const auto n = aMBS->q_.size();
		const auto indCoordsInd = dynSimul.independent_coordinate_indices();
		const auto d = dynSimul.independent_coordinate_indices().size();
		EXPECT_GT(d, 0);
		// const auto m = aMBS->Phi_q_.getNumRows();

		auto noise_dyn = gtsam::noiseModel::Isotropic::Sigma(d, 0.1);

		gtsam::Values valuesForQ;

		// Create a dummy factor:
		auto factorDyn = boost::make_shared<FactorDynamicsIndep>(
			&dynSimul, noise_dyn, Z(1), DZ(1), DDZ(1), Q(1), valuesForQ);

		// For different instants of time and mechanism positions and
		// velocities, test the factor jacobian:
		const double t_end = 3.0;  // end simulation time
		const double t_steps = 1.0;  // "large steps" to run the tests at

		dynSimul.params.time_step = 0.001;  // integrators timesteps

		mrpt::system::CTimeLogger timlog;
		timlog.enable(false);

		for (double t = 0; t < t_end;)
		{
			const double t_next = t + t_steps;
			dynSimul.run(t, t_next);
			t = t_next;

			std::cout << "Evaluating test for t=" << t << "\n";
			std::cout << "q  =" << aMBS->q_.transpose() << "\n";
			std::cout << "dq =" << aMBS->dotq_.transpose() << "\n";
			std::cout << "ddq =" << aMBS->ddotq_.transpose() << "\n";

			// Convert plain Eigen vectors into state_t classes (used as Values
			// in GTSAM factor graphs):
			const state_t q = state_t(aMBS->q_);
			const state_t z = state_t(mbse::subset(aMBS->q_, indCoordsInd));
			const state_t dotz =
				state_t(mbse::subset(aMBS->dotq_, indCoordsInd));
			const state_t ddotz =
				state_t(mbse::subset(aMBS->ddotq_, indCoordsInd));

			if (valuesForQ.exists(Q(1)))
				valuesForQ.at<mbse::state_t>(Q(1)) = q;
			else
				valuesForQ.insert(Q(1), q);

			// Evaluate theoretical Jacobians:
			gtsam::Matrix H[3];
			timlog.enter("factorsDynIndep.theoretical_jacob");

			factorDyn->evaluateError(z, dotz, ddotz, H[0], H[1], H[2]);

			timlog.leave("factorsDynIndep.theoretical_jacob");

			// Evaluate numerical Jacobians:
			gtsam::Matrix H_num[3];

			timlog.enter("factorsDyn.Indepnumeric_jacob");
			for (int i = 0; i < 3; i++)
			{
				NumericJacobParams p;
				p.q0 = q;
				p.z = z;
				p.dz = dotz;
				p.ddz = ddotz;
				p.diff_variable = i;
				p.factor = factorDyn.get();

				const gtsam::Vector& x = i == 0 ? p.z : (i == 1 ? p.dz : p.ddz);
				const gtsam::Vector x_incr =
					Eigen::VectorXd::Constant(x.rows(), x.cols(), 1e-9);

				mrpt::math::estimateJacobian(
					x,
					std::function<void(
						const gtsam::Vector& new_q, const NumericJacobParams& p,
						gtsam::Vector& err)>(&num_err_wrt_state),
					x_incr, p, H_num[i]);

				// Check:
				EXPECT_NEAR(
					(H[i] - H_num[i]).array().abs().maxCoeff(), 0.0, 1e-2)
					<< "H[" << i << "] Theoretical:\n"
					<< H[i]
					<< "\n"
					   "H_num["
					<< i << "] Numerical:\n"
					<< H_num[i] << "\n";
			}
			timlog.leave("factorsDyn.Indepnumeric_jacob");
		}
	}
	catch (const std::exception& e)
	{
		throw std::runtime_error(mrpt::exception_to_str(e));
	}
}
