/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// Example of dynamics
// ------------------------------------------------------------
#include <fstream>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/NonlinearEquality.h>
#include <gtsam_unstable/nonlinear/BatchFixedLagSmoother.h>
#include <gtsam_unstable/nonlinear/IncrementalFixedLagSmoother.h>
#include <iostream>
#include <mbse/AssembledRigidModel.h>
#include <mbse/ModelDefinition.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/factors/FactorConstraints.h>
#include <mbse/factors/FactorConstraintsVel.h>
#include <mbse/factors/FactorDynamics.h>
#include <mbse/factors/FactorEulerInt.h>
#include <mbse/factors/FactorTrapInt.h>
#include <mbse/model-examples.h>

#include <mrpt/3rdparty/tclap/CmdLine.h>
#include <mrpt/system/os.h>

// Declare the supported command line switches ===========
TCLAP::CmdLine cmd("smoother-forward-dynamics", ' ');

TCLAP::ValueArg<std::string> arg_mechanism(
	"", "mechanism", "Mechanism model YAML file", true, "mechanism.yaml",
	"YAML model definition", cmd);
TCLAP::ValueArg<double> arg_step_time(
	"", "dt", "Step time", false, 0.005, "xxx", cmd);
TCLAP::ValueArg<double> arg_end_time(
	"", "end-time", "Simulation end time", false, 5.0, "xxx", cmd);
TCLAP::ValueArg<double> arg_lag_time(
	"", "lag-time", "Smoother lag time", false, 0.1, "xxx", cmd);
TCLAP::ValueArg<unsigned int> arg_smoother_iterations(
	"", "smoother-iterations", "Smoother optimization maximum iterations",
	false, 15, "xxx", cmd);
TCLAP::ValueArg<std::string> arg_output_prefix(
	"", "output-prefix", "Output files prefix", false, "", "prefix", cmd);
TCLAP::SwitchArg arg_dont_add_q_constraints(
	"", "dont-add-q-constraints",
	"Do NOT add the q manifold constraint factors", cmd);
TCLAP::SwitchArg arg_dont_add_dq_constraints(
	"", "dont-add-dq-constraints",
	"Do NOT add the dq manifold constraint factors",

	cmd);
TCLAP::SwitchArg arg_show_factor_errors(
	"", "show-factor-errors", "Show factor errors for the final state", cmd);

TCLAP::SwitchArg argRunFinalBatch(
	"", "final-batch", "Run an additional final batch optimizer", cmd);

TCLAP::SwitchArg arg_verbose("v", "verbose", "Verbose console output", cmd);

void test_smoother()
{
	using gtsam::symbol_shorthand::A;
	using gtsam::symbol_shorthand::Q;
	using gtsam::symbol_shorthand::V;
	using namespace mbse;

	// Create the multibody object:
	const auto modelYaml =
		mrpt::containers::yaml::FromFile(arg_mechanism.getValue());
	const ModelDefinition model = ModelDefinition::FromYAML(modelYaml);

	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	aMBS->printCoordinates();

	CDynamicSimulator_R_matrix_dense dynSimul(aMBS);

	// Must be called before solve_ddotq(), needed inside the dynamics factors
	dynSimul.prepare();

	// iSAM2: Incrementally build the graph:
	// gtsam::NonlinearFactorGraph new_factors;
	//  gtsam::Values new_values;

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->q_.size();
	const auto m = aMBS->Phi_q_.getNumRows();
	std::cout << "problem n=" << n << " m=" << m << "\n";

	// x1, *y1*, x2, y2
	// 0   1     2   3
	std::vector<size_t> indepCoordIndices;
	indepCoordIndices.push_back(0);

	// Velocity prior: large sigma for all dq(i), except dq(i_indep)
	gtsam::Vector prior_dq_sigmas;
	const double large_std = 1;
	const double small_std = 1e-3;
	prior_dq_sigmas.setConstant(n, large_std);
	for (auto idx : indepCoordIndices) prior_dq_sigmas(idx) = small_std;

	auto noise_prior_dq_0 =
		gtsam::noiseModel::Diagonal::Sigmas(prior_dq_sigmas);

	auto noise_dyn = gtsam::noiseModel::Isotropic::Sigma(n, small_std);
	auto noise_constr_q = gtsam::noiseModel::Isotropic::Sigma(m, small_std);
	auto noise_constr_dq = gtsam::noiseModel::Isotropic::Sigma(m, small_std);
	auto noise_vel = gtsam::noiseModel::Isotropic::Sigma(n, small_std);
	auto noise_acc = gtsam::noiseModel::Isotropic::Sigma(n, small_std);

	const double dt = arg_step_time.getValue();
	const double t_end = arg_end_time.getValue();
	double t = 0;
	unsigned int N = static_cast<unsigned int>(t_end / dt);

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Create a feasible Q(0):
	aMBS->q_.setZero();
	aMBS->dotq_.setZero();
	aMBS->ddotq_.setZero();

	AssembledRigidModel::TComputeDependentParams cdp;	// default params
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
	state_t last_q = q_0, last_dq = zeros, last_ddq = zeros;

	std::multimap<double, gtsam::NonlinearFactor::shared_ptr> factorsByTime;
	std::multimap<double, gtsam::Key> keysByTime;

	// Create Prior factors:
	factorsByTime.emplace(
		0.0, boost::make_shared<gtsam::NonlinearEquality<state_t>>(Q(0), q_0));
	factorsByTime.emplace(
		0.0, boost::make_shared<gtsam::PriorFactor<state_t>>(
				 V(0), zeros, noise_prior_dq_0));

	const double lag = arg_lag_time.getValue();	 // seconds

	gtsam::NonlinearFactorGraph wholeFG;
	gtsam::Values wholeValues;

#define USE_BATCH_OPTIMIZER 1

#if USE_BATCH_OPTIMIZER
	gtsam::BatchFixedLagSmoother smoother(lag);

	smoother.params().maxIterations = arg_smoother_iterations.getValue();
#else
	gtsam::IncrementalFixedLagSmoother smoother(lag);

	// maxIterations = arg_smoother_iterations.getValue();
	// smoother.params().optimizationParams =  gnparams;
#endif

	gtsam::Values estimated;
	gtsam::FixedLagSmoother::KeyTimestampMap new_timestamps;

#if 0
  new_values.insert(Q(0), last_q);
  new_values.insert(V(0), last_dq);
  new_values.insert(A(0), last_ddq);
#endif

	wholeValues.insert(Q(0), last_q);
	wholeValues.insert(V(0), last_dq);
	wholeValues.insert(A(0), last_ddq);

	new_timestamps[Q(0)] = 0 * dt;
	new_timestamps[V(0)] = 0 * dt;
	new_timestamps[A(0)] = 0 * dt;

	// Save states to files:
	mrpt::math::CMatrixDouble Qs(N + 1, n + 1), dotQs(N + 1, n + 1),
		ddotQs(N + 1, n + 1);

	auto lambda_Values_toQ_DQ_DDQ = [dt, &Qs, &dotQs, &ddotQs, &last_q,
									 &last_dq,
									 &last_ddq](const gtsam::Values& values) {
		for (const auto& kv : values)
		{
			gtsam::Symbol s(kv.key);
			const auto step = s.index();
			const state_t val = values.at<state_t>(s);

			switch (s.chr())
			{
				case 'q':
					Qs.block(step, 1, 1, val.rows()) = val.transpose();
					Qs(step, 0) = dt * step;
					last_q = val;
					break;
				case 'v':
					dotQs.block(step, 1, 1, val.rows()) = val.transpose();
					dotQs(step, 0) = dt * step;
					last_dq = val;
					break;
				case 'a':
					ddotQs.block(step, 1, 1, val.rows()) = val.transpose();
					ddotQs(step, 0) = dt * step;
					last_ddq = val;
					break;
			};
		}
	};

	for (unsigned int timeStep = 0; timeStep < N; timeStep++, t += dt)
	{
		// Create Trapezoidal Integrator factors:
		factorsByTime.emplace(
			t, boost::make_shared<FactorTrapInt>(
				   dt, noise_vel, Q(timeStep), Q(timeStep + 1), V(timeStep),
				   V(timeStep + 1)));
		factorsByTime.emplace(
			t, boost::make_shared<FactorTrapInt>(
				   dt, noise_acc, V(timeStep), V(timeStep + 1), A(timeStep),
				   A(timeStep + 1)));

		// Create Dynamics factors:
		factorsByTime.emplace(
			t + dt, boost::make_shared<FactorDynamics>(
						&dynSimul, noise_dyn, Q(timeStep + 1), V(timeStep + 1),
						A(timeStep + 1)));
		if (timeStep == 0)
			factorsByTime.emplace(
				t, boost::make_shared<FactorDynamics>(
					   &dynSimul, noise_dyn, Q(timeStep), V(timeStep),
					   A(timeStep)));

		// Add dependent-coordinates constraint factor:
		if (!arg_dont_add_q_constraints.isSet())
			factorsByTime.emplace(
				t, boost::make_shared<FactorConstraints>(
					   aMBS, noise_constr_q, Q(timeStep)));

		if (!arg_dont_add_dq_constraints.isSet())
			factorsByTime.emplace(
				t, boost::make_shared<FactorConstraintsVel>(
					   aMBS, noise_constr_dq, Q(timeStep), V(timeStep)));

		// Create initial estimates:
		if (timeStep > 0) last_ddq = wholeValues.at<state_t>(A(timeStep - 1));

		wholeValues.insert(A(timeStep + 1), last_ddq);
		new_timestamps[A(timeStep + 1)] = (timeStep + 1) * dt;

		// Create initial estimates (so we can run the optimizer)
		wholeValues.insert(Q(timeStep + 1), last_q);
		new_timestamps[Q(timeStep + 1)] = (timeStep + 1) * dt;

		wholeValues.insert(V(timeStep + 1), last_dq);
		new_timestamps[V(timeStep + 1)] = (timeStep + 1) * dt;

		keysByTime.emplace(t, Q(timeStep));
		keysByTime.emplace(t, V(timeStep));
		keysByTime.emplace(t, A(timeStep));
		keysByTime.emplace(t, Q(timeStep + 1));
		keysByTime.emplace(t, V(timeStep + 1));
		keysByTime.emplace(t, A(timeStep + 1));

#if 0
    mbse::timelog().enter("smoother.update");
    const auto res = smoother.update(new_factors, new_values, new_timestamps);
    mbse::timelog().leave("smoother.update");

    std::cout << "n=" << timeStep << "/" << N
              << " keys: " << smoother.timestamps().size() << " rmse: "
              << std::sqrt(res.getError() / smoother.getFactors().size())
              << "\n";

    mbse::timelog().enter("smoother.calculateEstimate");
    estimated = smoother.calculateEstimate();
    mbse::timelog().leave("smoother.calculateEstimate");

    // reset containers:
    new_factors = gtsam::NonlinearFactorGraph();
    new_values.clear();
    new_timestamps.clear();
#else
		// Optimize a sliding window of the whole FG:
		gtsam::NonlinearFactorGraph fgWindow;
		for (auto timFactor : factorsByTime)
			if (timFactor.first >= (t - lag)) fgWindow += timFactor.second;

		gtsam::Values valuesWindow;
		uint64_t firstTimeIndexInWindow = std::numeric_limits<uint64_t>::max();
		for (auto timKey : keysByTime)
			if (timKey.first >= (t - lag))
			{
				const auto symbolIdx = gtsam::Symbol(timKey.second).index();
				mrpt::keep_min(firstTimeIndexInWindow, symbolIdx);

				if (!valuesWindow.exists(timKey.second))
					valuesWindow.insert(
						timKey.second, wholeValues.at(timKey.second));
			}

		// For any timestep > 0, add anchoring prior factors to trust the first
		// q,dq values in the window. Not neccesary for t==0 just before we
		// already have prior factors defined (permanent ones, for t=0).
		if (firstTimeIndexInWindow > 0)
		{
			const state_t q_init_win =
				wholeValues.at<state_t>(Q(firstTimeIndexInWindow));
			const state_t dq_init_win =
				wholeValues.at<state_t>(V(firstTimeIndexInWindow));

			// Create Prior factors:
			fgWindow.emplace_shared<gtsam::NonlinearEquality<state_t>>(
				Q(firstTimeIndexInWindow), q_init_win);
			fgWindow.emplace_shared<gtsam::NonlinearEquality<state_t>>(
				V(firstTimeIndexInWindow), dq_init_win);
		}

		gtsam::LevenbergMarquardtParams lp =
			gtsam::LevenbergMarquardtParams::LegacyDefaults();

		lp.maxIterations = arg_smoother_iterations.getValue();

#if 0
    lp.iterationHook = [&fgWindow](size_t iter, double errBef,
                                   double errAfter) {
      const auto N = fgWindow.size();
      std::cout << "LM iter #" << iter << " rmse: " << std::sqrt(errBef / N)
                << " -> " << std::sqrt(errAfter / N) << std::endl;
    };
#endif
		gtsam::LevenbergMarquardtOptimizer lm(fgWindow, valuesWindow, lp);
		const auto& estimated = lm.optimize();

#if 0
    std::cout << " === INIT:\n";
    fgWindow.printErrors(valuesWindow);
    valuesWindow.print("valuesWindow");
    std::cout << " === FINAL:\n";
    fgWindow.printErrors(estimated);
    estimated.print("estimated:");
#endif

		const double errorBeforeLM = fgWindow.error(valuesWindow);
		const double errorAfterLM = fgWindow.error(estimated);
		const auto numFactorsLM = fgWindow.size();

		std::cout << "n=" << timeStep << "/" << N
				  << " sliding window LevMarq."
					 " ErrorBefore = "
				  << errorBeforeLM
				  << " RMSE=" << std::sqrt(errorBeforeLM / numFactorsLM)
				  << " ErrorAfter  = " << errorAfterLM
				  << " RMSE=" << std::sqrt(errorAfterLM / numFactorsLM)
				  << " numFactors=" << numFactorsLM
				  << " iters:" << lm.iterations() << "\n";

#endif

		// save/update the last N values (older are "more refined"):
		for (const auto& kv : estimated)
		{
			// Update list of all values (for final batch):
			if (wholeValues.find(kv.key) == wholeValues.end())
				wholeValues.insert(kv.key, kv.value);
			else
				wholeValues.update(kv.key, estimated.at<state_t>(kv.key));
		}
		//    wholeValues.print("wholeValues:");

		// Update values in vectors for saving to disk:
		lambda_Values_toQ_DQ_DDQ(estimated);
	}

	// Build FG with all factors:
	for (auto timFactor : factorsByTime) wholeFG += timFactor.second;

#if 0
  const double errorAfter = smoother.getFactors().error(estimated);
  const auto numFactors = smoother.getFactors().size();

  std::cout << "ErrorAfter=" << errorAfter
            << " RMSE=" << std::sqrt(errorAfter / numFactors)
            << " numFactors=" << numFactors << "\n";
#endif

	auto lmbdPrintErr = [](const gtsam::Factor* /*factor*/,
						   double whitenedError, size_t /*index*/) -> bool {
		return true;  // whitenedError > 1e-3;
	};

	if (arg_show_factor_errors.isSet())
	{
		std::cout << "======== FACTOR ERRORS ============\n";
		// smoother.getFactors().
		wholeFG.printErrors(
			wholeValues /*estimated*/,
			"FACTOR ERRORS: ", gtsam::DefaultKeyFormatter, lmbdPrintErr);
	}

	// Additional batch step?
	if (argRunFinalBatch.isSet())
	{
		gtsam::LevenbergMarquardtParams lp =
			gtsam::LevenbergMarquardtParams::LegacyDefaults();

		lp.iterationHook = [&wholeFG](
							   size_t iter, double errBef, double errAfter) {
			const auto N = wholeFG.size();
			std::cout << "LM iter #" << iter
					  << " rmse: " << std::sqrt(errBef / N) << " -> "
					  << std::sqrt(errAfter / N) << std::endl;
		};

		lp.print("LevMarq parameters:");

		std::cout << "\n=== Running a batch optimization pass ===\n";
		gtsam::LevenbergMarquardtOptimizer lm(wholeFG, wholeValues, lp);
		const auto& lmValues = lm.optimize();

		const double errorBeforeLM = wholeFG.error(wholeValues);
		const double errorAfterLM = wholeFG.error(lmValues);
		const auto numFactorsLM = wholeFG.size();
		std::cout << "GLOBAL LEVMARQ:\n"
					 " ErrorBefore = "
				  << errorBeforeLM
				  << " RMSE=" << std::sqrt(errorBeforeLM / numFactorsLM)
				  << " numFactors=" << numFactorsLM << "\n"
				  << " ErrorAfter  = " << errorAfterLM
				  << " RMSE=" << std::sqrt(errorAfterLM / numFactorsLM)
				  << " numFactors=" << numFactorsLM << "\n";

		// Update values in vectors for saving to disk:
		lambda_Values_toQ_DQ_DDQ(lmValues);

		if (arg_show_factor_errors.isSet())
		{
			wholeFG.printErrors(
				lmValues, "LM FACTOR ERRORS: ", gtsam::DefaultKeyFormatter,
				lmbdPrintErr);
		}
	}

	const auto prefix = arg_output_prefix.getValue();

	std::cout << "Saving results to TXT files with prefix '" << prefix
			  << "'...\n";

	using namespace std::string_literals;  // "s"

	Qs.saveToTextFile(prefix + "q.txt"s);
	dotQs.saveToTextFile(prefix + "dq.txt"s);
	ddotQs.saveToTextFile(prefix + "ddq.txt"s);
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
