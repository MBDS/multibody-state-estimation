/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// clang-format off
/* Example invokation:

cd build
bin/mbse-fg-inverse-dynamics   \
  --mechanism ../config/mechanisms/fourbars1-with-rel-angle.yaml \
  --desired-trajectory ../config/trajectories/fourbars1-with-rel-angle-trajectory.txt \
  --imposed-coordinates "[ 4 ]" 

bin/mbse-fg-inverse-dynamics  \
  --mechanism ../config/mechanisms/pick-and-place-robot.yaml \
  --desired-trajectory motors.txt \
  --imposed-coordinates "[ 20 ; 21 ]" \
  --verbose \
  --lm-iterations 15

bin/mbse-fg-inverse-dynamics   \
  --mechanism ../config/mechanisms/fourbars1-with-rel-angle.yaml \
  --desired-trajectory ../config/trajectories/fourbars1-with-rel-angle-trajectory.txt \
  --enforce-coordinate-precision 10 \
  --enforce-force-sigma 1e3 \
  --imposed-coordinates "[ 4 ]" \
  --noise-vel-sigma 1e-1 \
  --noise-acc-sigma 1e-1 \
  --vel-constraints-sigma 1 \
  --dynamics-sigma 1 \
  --pos-constraints-sigma 1

 */
// clang-format on

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
	"", "desired-trajectory",
	"Desired mechanism trajectory, with a first column with timestamps, and "
	"one extra column per imposed degree-of-freedom. The number of these "
	"additional columns must match the number of coordinates given in "
	"--imposed-coordinates",
	true, "", "file.txt", cmd);

TCLAP::ValueArg<double> arg_enforcement_prec(
	"", "enforce-coordinate-precision",
	"Desired mechanism trajectory enforcement precision", false, 10, "10.0",
	cmd);

TCLAP::ValueArg<double> arg_force_enforcement_sigma(
	"", "enforce-force-sigma",
	"A priori factors for null forces in non-actuated Q coordinates (sigma)",
	false, 1e3, "1e3", cmd);

TCLAP::ValueArg<std::string> arg_indep_coords(
	"", "imposed-coordinates",
	"Comma-separated list of indices of coordinates whose motion is imposed "
	"via the trajectory file",
	true, "", "[ 4 ]", cmd);

TCLAP::ValueArg<double> arg_gravity(
	"", "gravity", "Gravity acceleration in Y direction", false, -9.81, "g",
	cmd);

TCLAP::ValueArg<double> arg_dynamics_sigma(
	"", "dynamics-sigma", "Sigma for the dynamics factors", false, 1.0, "1.0",
	cmd);

TCLAP::ValueArg<double> arg_q_constr_sigma(
	"", "pos-constraints-sigma", "Sigma for the position constratint factors",
	false, 1.0, "1.0", cmd);

TCLAP::ValueArg<double> arg_dq_constr_sigma(
	"", "vel-constraints-sigma", "Sigma for the velocity constratint factors",
	false, 1.0, "1.0", cmd);

TCLAP::ValueArg<double> arg_noise_vel_sigma(
	"", "noise-vel-sigma", "Sigma for q-dq integration", false, 0.1, "0.1",
	cmd);
TCLAP::ValueArg<double> arg_noise_acc_sigma(
	"", "noise-acc-sigma", "Sigma for dq-ddq integration", false, 0.1, "0.1",
	cmd);

TCLAP::ValueArg<unsigned int> arg_lm_iterations(
	"", "lm-iterations", "LevMarq optimization maximum iterations", false, 100,
	"xxx", cmd);

TCLAP::SwitchArg arg_verbose("v", "verbose", "Verbose console output", cmd);

TCLAP::SwitchArg arg_skipInverseDynamics(
	"", "skip-inverse-dynamics",
	"Run all preliminary steps but skip actual inverse dynamics, saving the "
	"preliminary results.",
	cmd);

TCLAP::SwitchArg arg_debugPrintGraphs(
	"", "debug-print-graphs", "Debug: print factor graphs", cmd);

TCLAP::ValueArg<double> arg_printFactorErrors(
	"", "print-factor-errors",
	"Print factor errors for those with error above the given constraints",
	false, 1.0, "1.0", cmd);

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

	const unsigned int N = trajectory.rows();
	std::cout << "Desired trajectory loaded with " << N << " time points.\n";

	// autodetect timestep:
	const double dt = trajectory(1, 0) - trajectory(0, 0);
	std::cout << "Timestep (from trajectory file): " << dt << std::endl;
	trajectory.removeColumns({0});

	auto noise_vel =
		gtsam::noiseModel::Isotropic::Sigma(n, arg_noise_vel_sigma.getValue());
	auto noise_acc =
		gtsam::noiseModel::Isotropic::Sigma(n, arg_noise_acc_sigma.getValue());

	// Indep coords:
	mrpt::math::CMatrixDouble indepCoordsVec;
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

	const size_t nImposedDOFs = indepCoordIndices.size();

	ASSERT_GE_(trajectory.cols(), indepCoordsVec.size());

	// Enforcement noise model:
	mrpt::math::CVectorDouble posEnforcePrecisions;
	posEnforcePrecisions.setZero(n);
	for (const auto i : indepCoordIndices)
		posEnforcePrecisions[i] = arg_enforcement_prec.getValue();
	const gtsam::Vector posEnforcePrecisionsEig =
		posEnforcePrecisions.asEigen();
	const auto noise_pos_enforcement =
		gtsam::noiseModel::Diagonal::Precisions(posEnforcePrecisionsEig);

	// Enforcement noise model for forces:
	mrpt::math::CVectorDouble forceZerosPrioriEnforcementSigmas;
	forceZerosPrioriEnforcementSigmas.setZero(n);
	for (const auto i : indepCoordIndices)
		forceZerosPrioriEnforcementSigmas[i] =
			arg_force_enforcement_sigma.getValue();
	const gtsam::Vector forceZerosPrioriEnforcementSigmasEig =
		forceZerosPrioriEnforcementSigmas.asEigen();
	const auto forceZerosPrioriEnforcement =
		gtsam::noiseModel::Diagonal::Sigmas(
			forceZerosPrioriEnforcementSigmasEig);

	// Between Q factor for smooth motion (solve branching).
	gtsam::Vector betweenQSigmas;
	const double large_std = 1;
	betweenQSigmas.setConstant(n, large_std);

	auto noise_between_q = gtsam::noiseModel::Diagonal::Sigmas(betweenQSigmas);
	auto noise_dyn =
		gtsam::noiseModel::Isotropic::Sigma(n, arg_dynamics_sigma.getValue());
	auto noise_constr_q =
		gtsam::noiseModel::Isotropic::Sigma(m, arg_q_constr_sigma.getValue());
	auto noise_constr_dq =
		gtsam::noiseModel::Isotropic::Sigma(m, arg_dq_constr_sigma.getValue());
	auto noise_constant_F = gtsam::noiseModel::Isotropic::Sigma(n, 1.0);

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Extract an initial q_0 from the assembled multibody problem:
	std::cout << "nImposedDOFs: " << nImposedDOFs << "\n";

	// Position enforcement factor:
	{
		gtsam::Vector qn = gtsam::Vector::Zero(n);
		for (size_t i = 0; i < nImposedDOFs; i++)
			qn[indepCoordIndices.at(i)] = trajectory(0, i);

		aMBS->q_ = qn;
		AssembledRigidModel::ComputeDependentParams dp;
		dp.maxPhiNorm = 1e-8;
		dp.nItersMax = 100;
		AssembledRigidModel::ComputeDependentResults dr;

		aMBS->computeDependentPosVelAcc(indepCoordIndices, true, false, dp, dr);

		std::cout << "Position problem final |Phi(q)|=" << aMBS->Phi_.norm()
				  << ", Phi(q)=" << aMBS->Phi_.transpose() << std::endl;

		ASSERT_LT_(dr.pos_final_phi, 0.01);
	}

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
		gtsam::Vector qn = gtsam::Vector::Zero(n);
		for (size_t i = 0; i < nImposedDOFs; i++)
			qn[indepCoordIndices.at(i)] = trajectory(timeStep, i);

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
				Q(timeStep - 1), Q(timeStep), zeros, noise_between_q);
		}
	}

	gtsam::LevenbergMarquardtParams optParams;
	// optParams.verbosityLM = gtsam::LevenbergMarquardtParams::DAMPED;
	optParams.lambdaUpperBound = 1e10;
	optParams.lambdaFactor = 2.0;
	optParams.diagonalDamping = true;
	optParams.absoluteErrorTol = 0;
	optParams.relativeErrorTol = 1e-5;
	optParams.maxIterations = arg_lm_iterations.getValue();

	if (arg_verbose.isSet())
	{
		optParams.iterationHook = [](size_t iter, double oldError,
									 double newError) {
			std::cout << "- LM iteration " << iter << " error: " << oldError
					  << " -> " << newError << std::endl;
		};
	}

	const auto lmbdPrintErr = [](const gtsam::Factor*, double err, size_t) {
		return err > arg_printFactorErrors.getValue();
	};
	const auto defKeyFrm = gtsam::DefaultKeyFormatter;

	{
		const auto numFactors = fg.size();
		const double errorBefore = fg.error(values);

		std::cout << " PASS 1: q only\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, optParams);
		if (arg_debugPrintGraphs.isSet()) fg.print();

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";
		values = result1;

		if (arg_printFactorErrors.isSet())
			fg.printErrors(values, "", defKeyFrm, lmbdPrintErr);
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

		std::cout << " PASS 2: q,dq,ddq\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, optParams);

		if (arg_debugPrintGraphs.isSet()) fg.print();

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;

		if (arg_printFactorErrors.isSet())
			fg.printErrors(values, "", defKeyFrm, lmbdPrintErr);
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

		std::cout << "  PASS 3: velocity constraints\n"
					 " ==================================\n";

		std::cout << " ErrorBefore=" << errorBefore
				  << " RMSE=" << std::sqrt(errorBefore / numFactors)
				  << " numFactors=" << numFactors << "\n";

		// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

		gtsam::LevenbergMarquardtOptimizer lm1(fg, values, optParams);

		if (arg_debugPrintGraphs.isSet()) fg.print();

		const gtsam::Values& result1 = lm1.optimize();

		const double errorAfter = fg.error(result1);

		std::cout << " ErrorAfter=" << errorAfter
				  << " RMSE=" << std::sqrt(errorAfter / numFactors)
				  << " numFactors=" << numFactors
				  << " Iterations: " << lm1.iterations() << "\n\n";

		values = result1;

		if (arg_printFactorErrors.isSet())
			fg.printErrors(values, "", defKeyFrm, lmbdPrintErr);
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
	if (!arg_skipInverseDynamics.isSet())
	{
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

			std::cout << " PASS 4: inverse dynamics\n"
						 " ==================================\n";

			std::cout << " ErrorBefore=" << errorBefore
					  << " RMSE=" << std::sqrt(errorBefore / numFactors)
					  << " numFactors=" << numFactors << "\n";

			// lmParams.verbosityLM = gtsam::LevenbergMarquardtParams::LAMBDA;

			gtsam::LevenbergMarquardtOptimizer lm1(fg, values, optParams);

			if (arg_debugPrintGraphs.isSet()) fg.print();

			const gtsam::Values& result1 = lm1.optimize();

			const double errorAfter = fg.error(result1);

			std::cout << " ErrorAfter=" << errorAfter
					  << " RMSE=" << std::sqrt(errorAfter / numFactors)
					  << " numFactors=" << numFactors
					  << " Iterations: " << lm1.iterations() << "\n\n";

			values = result1;

			if (arg_printFactorErrors.isSet())
				fg.printErrors(values, "", defKeyFrm, lmbdPrintErr);
		}
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
	Qs.saveToTextFile("q.txt", {}, false, "% TIMESTAMP q[0]  ... q[n]\n");
	dotQs.saveToTextFile("dq.txt", {}, false, "% TIMESTAMP dq[0]  ... dq[n]\n");
	ddotQs.saveToTextFile(
		"ddq.txt", {}, false, "% TIMESTAMP ddq[0]  ... ddq[n]\n");
	Fs.saveToTextFile(
		"Q_forces.txt", {}, false, "% TIMESTAMP Q[0]  ... Q[n]\n");
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
