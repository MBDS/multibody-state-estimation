/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

/* Example invokation:

./mbse-fg-inverse-dynamics --desired-trajectory trajectory.txt    \
  --mechanism 4bars --enforce-coordinate-precisions "[0; 0; 0; 0; 1e3]"\
  --enforce-force-sigmas "[0; 0; 0; 0; 1e3]" \
  --indep-coord-indices "[ 4 ]" \
  --lm-iterations 300 --dt 5e-3

 */

#include <fstream>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>
#include <iostream>
#include <mbse/AssembledRigidModel.h>
#include <mbse/ModelDefinition.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/factors/FactorConstraints.h>
#include <mbse/factors/FactorConstraintsVel.h>
#include <mbse/factors/FactorEulerInt.h>
#include <mbse/factors/FactorInverseDynamics.h>
#include <mbse/factors/FactorTrapInt.h>
#include <mbse/model-examples.h>
#include <mrpt/core/round.h>
#include <mrpt/math/CVectorDynamic.h>

#include <mrpt/3rdparty/tclap/CmdLine.h>

// Declare the supported command line switches ===========
TCLAP::CmdLine cmd("mbse-fg-inverse-dynamics", ' ');

TCLAP::ValueArg<std::string> arg_mechanism(
	"", "mechanism", "Mechanism model YAML file", true, "mechanism.yaml",
	"YAML model definition", cmd);
TCLAP::ValueArg<std::string> arg_desired_trajectory(
	"", "desired-trajectory", "Desired mechanism trajectory", true, "",
	"file.txt", cmd);
TCLAP::ValueArg<std::string> arg_enforcement_prec(
	"", "enforce-coordinate-precisions",
	"Desired mechanism trajectory enforcement precisions", true, "",
	"[1e3 ... 0 0]", cmd);

TCLAP::ValueArg<std::string> arg_force_enforcement_sigmas(
	"", "enforce-force-sigmas", "A priori factors for forces in Q (sigmas)",
	true, "", "[1e3 ... 0 0]", cmd);

TCLAP::ValueArg<std::string> arg_indep_coords(
	"", "indep-coord-indices", "Indices of independent coordinates", true, "",
	"[ 4 ]", cmd);

TCLAP::ValueArg<double> arg_step_time(
	"", "dt", "Step time", false, 0.005, "xxx", cmd);
TCLAP::ValueArg<double> arg_gravity(
	"", "gravity", "Gravity acceleration in Y direction", false, -9.81, "g",
	cmd);
TCLAP::ValueArg<unsigned int> arg_lm_iterations(
	"", "lm-iterations", "LevMarq optimization maximum iterations", false, 100,
	"xxx", cmd);

TCLAP::SwitchArg arg_verbose("v", "verbose", "Verbose console output", cmd);

void test_smoother()
{
	using gtsam::symbol_shorthand::A;  // ddq_k
	using gtsam::symbol_shorthand::F;  // Q_k
	using gtsam::symbol_shorthand::Q;  // q_k
	using gtsam::symbol_shorthand::V;  // dq_k
	using namespace mbse;

	// Create the multibody object:
	const auto modelYaml =
		mrpt::containers::yaml::FromFile(arg_mechanism.getValue());
	const ModelDefinition model = ModelDefinition::FromYAML(modelYaml);

	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();

	aMBS->setGravityVector(0, arg_gravity.getValue(), 0);

	aMBS->printCoordinates();

	CDynamicSimulator_R_matrix_dense dynSimul(aMBS);
	// CDynamicSimulator_Lagrange_LU_dense dynSimul(aMBS);

	// Must be called before solve_ddotq(), needed inside the dynamics factors
	dynSimul.prepare();

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->q_.size();
	const auto m = aMBS->Phi_q_.getNumRows();
	std::cout << "problem n=" << n << " m=" << m << "\n";

	// Load enforced trajectory:
	mrpt::math::CMatrixDouble trajectory;
	trajectory.loadFromTextFile(arg_desired_trajectory.getValue());
	ASSERT_EQUAL_(trajectory.cols(), n);

	const unsigned int N = trajectory.rows();
	std::cout << "enforcement trajectory loaded with " << N
			  << " time points.\n";

	const double noise_vel_sigma = 0.01, noise_acc_sigma = 0.01;

	auto noise_vel = gtsam::noiseModel::Isotropic::Sigma(n, noise_vel_sigma);
	auto noise_acc = gtsam::noiseModel::Isotropic::Sigma(n, noise_acc_sigma);

	// Enforcement noise model:
	mrpt::math::CVectorDouble posEnforcePrecisions;
	if (!posEnforcePrecisions.fromMatlabStringFormat(
			arg_enforcement_prec.getValue(), std::cerr))
		THROW_EXCEPTION_FMT(
			"Invalid matlab-like vector: '%s'",
			arg_enforcement_prec.getValue().c_str());
	std::cout << "q enforcement precisions: "
			  << posEnforcePrecisions.transpose() << "\n";
	const gtsam::Vector posEnforcePrecisionsEig =
		posEnforcePrecisions.asEigen();
	const auto noise_pos_enforcement =
		gtsam::noiseModel::Diagonal::Precisions(posEnforcePrecisionsEig);

	// Enforcement noise model for forces:
	mrpt::math::CVectorDouble forceZerosPrioriEnforcementSigmas;
	if (!forceZerosPrioriEnforcementSigmas.fromMatlabStringFormat(
			arg_force_enforcement_sigmas.getValue(), std::cerr))
		THROW_EXCEPTION_FMT(
			"Invalid matlab-like vector: '%s'",
			arg_force_enforcement_sigmas.getValue().c_str());
	std::cout << "Q enforcement sigmas: "
			  << forceZerosPrioriEnforcementSigmas.transpose() << "\n";
	const gtsam::Vector forceZerosPrioriEnforcementSigmasEig =
		forceZerosPrioriEnforcementSigmas.asEigen();
	const auto forceZerosPrioriEnforcement =
		gtsam::noiseModel::Diagonal::Sigmas(
			forceZerosPrioriEnforcementSigmasEig);

	// x1, y1, x2, y2,  theta
	// 0   1     2   3,  4
	mrpt::math::CVectorDouble indepCoordsVec;
	if (!indepCoordsVec.fromMatlabStringFormat(
			arg_indep_coords.getValue(), std::cerr))
		THROW_EXCEPTION_FMT(
			"Invalid matlab-like vector: '%s'",
			arg_indep_coords.getValue().c_str());
	std::cout << "Independent coordinate indices: "
			  << indepCoordsVec.transpose() << "\n";
	std::vector<size_t> indepCoordIndices;
	for (int i = 0; i < indepCoordsVec.size(); i++)
		indepCoordIndices.push_back(mrpt::round(indepCoordsVec[i]));

	// Velocity prior: large sigma for all dq(i), except dq(i_indep)
	gtsam::Vector prior_dq_sigmas, prior_ddq_sigmas;
	const double large_std = 1e3;
	prior_dq_sigmas.setConstant(n, large_std);
	prior_ddq_sigmas.setConstant(n, large_std);

	auto noise_prior_dq = gtsam::noiseModel::Diagonal::Sigmas(prior_dq_sigmas);
	auto noise_prior_ddq =
		gtsam::noiseModel::Diagonal::Sigmas(prior_ddq_sigmas);
	auto noise_dyn = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);
	auto noise_constr_q = gtsam::noiseModel::Isotropic::Sigma(m, 1e-3);
	auto noise_constr_dq = gtsam::noiseModel::Isotropic::Sigma(m, 1e-3);
	auto noise_constant_F = gtsam::noiseModel::Isotropic::Sigma(n, 1.0);

	const double dt = arg_step_time.getValue();

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Create a feasible Q(0):
	aMBS->q_.setZero();
	aMBS->dotq_.setZero();
	aMBS->ddotq_.setZero();

	AssembledRigidModel::TComputeDependentParams cdp;  // default params
	AssembledRigidModel::TComputeDependentResults cdr;
	// Solve the position problem:
	aMBS->q_[0] = 1;
	aMBS->q_[1] = 0.1;
	aMBS->q_[3] = 5;

	aMBS->computeDependentPosVelAcc(indepCoordIndices, true, true, cdp, cdr);
	std::cout << "Position problem final |Phi(q)|=" << cdr.pos_final_phi
			  << "\n";
	ASSERT_BELOW_(cdr.pos_final_phi, 1e-4);

	// Extract q_ from the assembled multibody problem:
	state_t q_0 = gtsam::Vector(aMBS->q_);
	std::cout << "q0: " << q_0.transpose() << "\n";

	/* =======================================================================
	 * PASS 1
	 *
	 * - Variables: q only
	 *
	 * - q priors with desired trajectory
	 * - position constraints
	 * - between factors for smooth motion between consecutive timesteps
	 *
	 * =======================================================================
	 */
	gtsam::NonlinearFactorGraph fg;
	gtsam::Values values;

	for (unsigned int timeStep = 0; timeStep < N; timeStep++)
	{
		// Position enforcement factor:
		gtsam::Vector qn;
		trajectory.extractRow(timeStep, qn);
		fg.emplace_shared<gtsam::PriorFactor<state_t>>(
			Q(timeStep), qn, noise_pos_enforcement);

		// Initial value:
		values.insert(Q(timeStep), q_0);

		// Add dependent-coordinates constraint factor:
		fg.emplace_shared<FactorConstraints>(aMBS, noise_constr_q, Q(timeStep));

		// between factor: required to solve branch indeterminatiosn (e.g. a
		// 2-bar mechanism with 2 possible branches)
		if (timeStep > 0)
		{
			fg.emplace_shared<gtsam::BetweenFactor<state_t>>(
				Q(timeStep - 1), Q(timeStep), zeros, noise_prior_dq);
		}
	}

	{
		const auto numFactors = fg.size();
		const double errorBefore = fg.error(values);

		std::cout << "             PASS 1\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		gtsam::LevenbergMarquardtParams lmParams;
		lmParams.maxIterations = arg_lm_iterations.getValue();
		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, lmParams);

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;
	}

	/* =======================================================================
	 * PASS 2
	 *
	 * - Variables: q, dq, ddq
	 *
	 * Extra factors:
	 * - time integration: q, dq
	 * - time integration: dq, ddq
	 *
	 * =======================================================================
	 */
	for (unsigned int timeStep = 0; timeStep < N; timeStep++)
	{
		// Create Trapezoidal Integrator factors:
		if (timeStep < N - 1)
		{
			fg.emplace_shared<FactorTrapInt>(
				dt, noise_vel, Q(timeStep), Q(timeStep + 1), V(timeStep),
				V(timeStep + 1));
			fg.emplace_shared<FactorTrapInt>(
				dt, noise_acc, V(timeStep), V(timeStep + 1), A(timeStep),
				A(timeStep + 1));
		}

		// Create initial estimates (so we can run the optimizer)
		values.insert(A(timeStep), zeros);
		values.insert(V(timeStep), zeros);
	}

	{
		const auto numFactors = fg.size();
		const double errorBefore = fg.error(values);

		std::cout << "             PASS 2\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		gtsam::LevenbergMarquardtParams lmParams;
		lmParams.maxIterations = arg_lm_iterations.getValue();
		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, lmParams);

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;
	}

	/* =======================================================================
	 * PASS 3
	 *
	 * - Variables: q, dq, ddq
	 *
	 * Extra factors:
	 * - velocity constraints
	 *
	 * =======================================================================
	 */
	for (unsigned int timeStep = 0; timeStep < N; timeStep++)
	{
		fg.emplace_shared<FactorConstraintsVel>(
			aMBS, noise_constr_dq, Q(timeStep), V(timeStep));
	}

	{
		const auto numFactors = fg.size();
		const double errorBefore = fg.error(values);

		std::cout << "             PASS 3\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		gtsam::LevenbergMarquardtParams lmParams;
		lmParams.maxIterations = arg_lm_iterations.getValue();
		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, lmParams);

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;
	}

	/* =======================================================================
	 * PASS 4
	 *
	 * - Variables: q, dq, ddq, Q (force)
	 *
	 * Extra factors:
	 * - dynamics factor
	 *
	 * =======================================================================
	 */

	for (unsigned int timeStep = 0; timeStep < N; timeStep++)
	{
		// A priori factor for forces:
		fg.emplace_shared<gtsam::PriorFactor<state_t>>(
			F(timeStep), zeros, forceZerosPrioriEnforcement);

		// Create Inverse Dynamics factors:
		fg.emplace_shared<FactorInverseDynamics>(
			&dynSimul, noise_dyn, Q(timeStep), V(timeStep), A(timeStep),
			F(timeStep));

		// Create initial estimates (so we can run the optimizer)
		values.insert(F(timeStep), zeros);
	}

	{
		const auto numFactors = fg.size();
		const double errorBefore = fg.error(values);

		std::cout << "             PASS 4\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		gtsam::LevenbergMarquardtParams lmParams;
		lmParams.maxIterations = arg_lm_iterations.getValue();
		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, lmParams);

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;
	}

	/* =======================================================================
	 * Extra values to matrices
	 * =======================================================================
	 */
	// Save states to files:
	mrpt::math::CMatrixDouble Qs(N, n + 1), dotQs(N, n + 1), ddotQs(N, n + 1),
		Fs(N, n + 1);

	for (const auto& kv : values)
	{
		gtsam::Symbol s(kv.key);
		const auto step = s.index();
		const state_t val = kv.value.cast<mbse::state_t>();
		switch (s.chr())
		{
			case 'q':
				Qs.block(step, 1, 1, val.rows()) = val.transpose();
				Qs(step, 0) = dt * step;
				break;
			case 'v':
				dotQs.block(step, 1, 1, val.rows()) = val.transpose();
				dotQs(step, 0) = dt * step;
				break;
			case 'a':
				ddotQs.block(step, 1, 1, val.rows()) = val.transpose();
				ddotQs(step, 0) = dt * step;
				break;
			case 'f':
				Fs.block(step, 1, 1, val.rows()) = val.transpose();
				Fs(step, 0) = dt * step;
				break;
		};
	}

	std::cout << "Saving results to TXT files...\n";
	Qs.saveToTextFile("q.txt");
	dotQs.saveToTextFile("dq.txt");
	ddotQs.saveToTextFile("ddq.txt");
	Fs.saveToTextFile("Q_forces.txt");
}

int main(int argc, char** argv)
{
	try
	{
		// Parse arguments:
		if (!cmd.parse(argc, argv))
			throw std::runtime_error("");  // should exit.

		test_smoother();
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << mrpt::exception_to_str(e);
	}
}
