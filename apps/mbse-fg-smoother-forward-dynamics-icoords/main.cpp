/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

/*
 Example of usage:

 ./smoother-forward-dynamics-icoord  --lag-time 5e-3 --dt 1e-3 \
	--mechanism  4bars \
	--indep-coord-indices "[ 4 ]"

*/

// Example of dynamics
// ------------------------------------------------------------
#include <fstream>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam_unstable/nonlinear/BatchFixedLagSmoother.h>
#include <gtsam_unstable/nonlinear/IncrementalFixedLagSmoother.h>
#include <iostream>
#include <mbse/AssembledRigidModel.h>
#include <mbse/ModelDefinition.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/factors/FactorConstraintsAccIndep.h>
#include <mbse/factors/FactorConstraintsIndep.h>
#include <mbse/factors/FactorConstraintsVelIndep.h>
#include <mbse/factors/FactorDynamicsIndep.h>
#include <mbse/factors/FactorEulerInt.h>
#include <mbse/factors/FactorTrapInt.h>
#include <mbse/mbse-utils.h>
#include <mrpt/core/round.h>
#include <mrpt/math/CVectorDynamic.h>

#include <mrpt/3rdparty/tclap/CmdLine.h>
#include <mrpt/system/os.h>

// Declare the supported command line switches ===========
TCLAP::CmdLine cmd("smoother-forward-dynamics-icoord", ' ');

TCLAP::ValueArg<std::string> arg_mechanism(
	"", "mechanism", "Mechanism model", true, "", "xxx", cmd);
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
TCLAP::SwitchArg arg_show_factor_errors(
	"", "show-factor-errors", "Show factor errors for the final state", cmd);

TCLAP::SwitchArg argRunFinalBatch(
	"", "final-batch", "Run an additional final batch optimizer", cmd);

TCLAP::SwitchArg arg_do_not_show_error_progress(
	"", "dont-show-error-progress",
	"Do NOT show the error for each timestep (this saves some CPU time).", cmd);

TCLAP::ValueArg<std::string> arg_indep_coords(
	"", "indep-coord-indices", "Indices of independent coordinates", true, "",
	"[ 4 ]", cmd);

TCLAP::SwitchArg arg_verbose("v", "verbose", "Verbose console output", cmd);

// Symbols:
constexpr auto sQ = gtsam::SymbolGenerator('q');
constexpr auto sQp = gtsam::SymbolGenerator('v');
constexpr auto sQpp = gtsam::SymbolGenerator('a');
constexpr auto sZ = gtsam::SymbolGenerator('n');
constexpr auto sZp = gtsam::SymbolGenerator('m');
constexpr auto sZpp = gtsam::SymbolGenerator('l');

void test_smoother()
{
	using namespace mbse;

	// Create the multibody object:
	const auto modelYaml =
		mrpt::containers::yaml::FromFile(arg_mechanism.getValue());
	const ModelDefinition model = ModelDefinition::FromYAML(modelYaml);

	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_Indep_dense dynSimul(aMBS);

	// x1, y1, x2, y2, theta
	// 0   1     2   3, 4
	// Enforce the use of the theta angle as independent coordinate:
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

	dynSimul.independent_coordinate_indices(indepCoordIndices);

	// Must be called before solve_ddotq(), needed inside the dynamics factors
	dynSimul.prepare();

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->q_.size();
	const auto m = aMBS->Phi_q_.getNumRows();
	std::cout << "problem n=" << n << " m=" << m << "\n";

	const auto nDofs = indepCoordIndices.size();
	std::cout << "problem: size(q)=" << n << " size(Phi)=" << m
			  << " size(z)=" << nDofs << "\n";
	aMBS->printCoordinates();

	auto noise_prior_z_0 =
		gtsam::noiseModel::Isotropic::Sigma(nDofs, 1e-10);	// fixed
	auto noise_prior_dz_0 =
		gtsam::noiseModel::Isotropic::Sigma(nDofs, 1e-10);	// fixed

	const double small_std = 1e-3;

	auto noise_dyn_z = gtsam::noiseModel::Isotropic::Sigma(nDofs, small_std);
	auto noise_constr_z =
		gtsam::noiseModel::Isotropic::Sigma(m + nDofs, small_std);
	auto noise_constr_dz =
		gtsam::noiseModel::Isotropic::Sigma(m + nDofs, small_std);
	auto noise_vel_z = gtsam::noiseModel::Isotropic::Sigma(nDofs, small_std);
	auto noise_acc_z = gtsam::noiseModel::Isotropic::Sigma(nDofs, small_std);
	const auto softBetweenNoise = gtsam::noiseModel::Isotropic::Sigma(n, 1e+2);

	const double dt = arg_step_time.getValue();
	const double t_end = arg_end_time.getValue();
	double t = 0;
	unsigned int N = static_cast<unsigned int>(t_end / dt);

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros_q = gtsam::Vector(gtsam::Vector::Zero(n, 1));
	const state_t zeros_z = gtsam::Vector(gtsam::Vector::Zero(nDofs, 1));

	// Create a feasible Q(0):
	if (0)
	{
		aMBS->q_.setZero();
		AssembledRigidModel::ComputeDependentParams cdp;  // default params
		AssembledRigidModel::ComputeDependentResults cdr;
		aMBS->computeDependentPosVelAcc(
			indepCoordIndices, true, true, cdp, cdr);
		std::cout << "Position problem final |Phi(q)|=" << cdr.pos_final_phi
				  << "\n";
		ASSERT_LT_(cdr.pos_final_phi, 1e-4);
	}
	else
	{
		aMBS->dotq_.setZero();
		aMBS->ddotq_.setZero();
	}

	// Extract q_ from the assembled multibody problem:
	const state_t q_0 = gtsam::Vector(aMBS->q_);
	const state_t z_0 = mbse::subset(q_0, indepCoordIndices);

	std::cout << "q0: " << q_0.transpose() << "\n";

	state_t last_q = q_0, last_dq = zeros_q, last_ddq = zeros_q;
	state_t last_z = z_0, last_dz = zeros_z, last_ddz = zeros_z;

	std::multimap<double, gtsam::NonlinearFactor::shared_ptr> factorsByTime;
	std::multimap<double, gtsam::Key> keysByTime;

	// Create Prior factors:
	factorsByTime.emplace(
		0.0, boost::make_shared<gtsam::PriorFactor<state_t>>(
				 sZ(0), z_0, noise_prior_z_0));
	factorsByTime.emplace(
		0.0, boost::make_shared<gtsam::PriorFactor<state_t>>(
				 sZp(0), zeros_z, noise_prior_dz_0));

	const double lag = arg_lag_time.getValue();	 // seconds

	gtsam::NonlinearFactorGraph wholeFG;
	gtsam::Values wholeValues;

	gtsam::Values estimated;
	gtsam::FixedLagSmoother::KeyTimestampMap new_timestamps;

	wholeValues.insert(sQ(0), last_q);
	new_timestamps[sQ(0)] = 0 * dt;

	wholeValues.insert(sQp(0), last_dq);
	new_timestamps[sQp(0)] = 0 * dt;

	wholeValues.insert(sQpp(0), last_ddq);
	new_timestamps[sQpp(0)] = 0 * dt;

	wholeValues.insert(sZ(0), last_z);
	new_timestamps[sZ(0)] = 0 * dt;

	wholeValues.insert(sZp(0), last_dz);
	new_timestamps[sZp(0)] = 0 * dt;

	wholeValues.insert(sZpp(0), last_ddz);
	new_timestamps[sZpp(0)] = 0 * dt;

	// Save states to files:
	mrpt::math::CMatrixDouble Qs(N + 1, n + 1), dotQs(N + 1, n + 1),
		ddotQs(N + 1, n + 1);

	auto lambda_Values_toQ_DQ_DDQ = [dt, &Qs, &dotQs, &ddotQs, &last_q,
									 &last_dq, &last_ddq, &last_z, &last_dz,
									 &last_ddz](const gtsam::Values& values) {
		for (const auto& kv : values)
		{
			gtsam::Symbol s(kv.key);
			const auto step = s.index();
			const state_t val = values.at<state_t>(s);

			switch (s.chr())
			{
				case sQ.chr():
					Qs.block(step, 1, 1, val.rows()) = val.transpose();
					Qs(step, 0) = dt * step;
					last_q = val;
					break;
				case sQp.chr():
					dotQs.block(step, 1, 1, val.rows()) = val.transpose();
					dotQs(step, 0) = dt * step;
					last_dq = val;
					break;
				case sQpp.chr():
					ddotQs.block(step, 1, 1, val.rows()) = val.transpose();
					ddotQs(step, 0) = dt * step;
					last_ddq = val;
					break;
				case sZ.chr():
					last_z = val;
					break;
				case sZp.chr():
					last_dz = val;
					break;
				case sZpp.chr():
					last_ddz = val;
					break;
			};
		}
	};

	const bool buildWholeFG =
		arg_show_factor_errors.isSet() || argRunFinalBatch.isSet();

	gtsam::LevenbergMarquardtParams lp =
		gtsam::LevenbergMarquardtParams::LegacyDefaults();

	lp.maxIterations = arg_smoother_iterations.getValue();
	lp.absoluteErrorTol = 0;
	lp.relativeErrorTol = 1e-8;

	for (unsigned int timeStep = 0; timeStep < N; timeStep++, t += dt)
	{
		mrpt::system::CTimeLoggerEntry tleStep(
			mbse::timelog(), "wholeTimeStep");

		// Create Trapezoidal Integrator factors:
		factorsByTime.emplace(
			t, boost::make_shared<FactorTrapInt>(
				   dt, noise_vel_z, sZ(timeStep), sZ(timeStep + 1),
				   sZp(timeStep), sZp(timeStep + 1)));
		factorsByTime.emplace(
			t, boost::make_shared<FactorTrapInt>(
				   dt, noise_acc_z, sZp(timeStep), sZp(timeStep + 1),
				   sZpp(timeStep), sZpp(timeStep + 1)));

		// Create Dynamics factors:
		factorsByTime.emplace(
			t + dt,
			boost::make_shared<FactorDynamicsIndep>(
				&dynSimul, noise_dyn_z, sZ(timeStep + 1), sZp(timeStep + 1),
				sZpp(timeStep + 1), sQ(timeStep + 1), wholeValues));
		if (timeStep == 0)
			factorsByTime.emplace(
				t, boost::make_shared<FactorDynamicsIndep>(
					   &dynSimul, noise_dyn_z, sZ(timeStep), sZp(timeStep),
					   sZpp(timeStep), sQ(timeStep), wholeValues));

		// "Soft equality" constraints between q_i and q_{i+1} to solve
		// configuration/branches ambiguities:
		factorsByTime.emplace(
			t, boost::make_shared<gtsam::BetweenFactor<state_t>>(
				   sQ(timeStep), sQ(timeStep + 1), zeros_q, softBetweenNoise));

		// Add dependent-coordinates constraint factor:
		if (timeStep == 0)
		{
			factorsByTime.emplace(
				t, boost::make_shared<FactorConstraintsIndep>(
					   aMBS, indepCoordIndices, noise_constr_z, sZ(timeStep),
					   sQ(timeStep)));

			factorsByTime.emplace(
				t, boost::make_shared<FactorConstraintsVelIndep>(
					   aMBS, indepCoordIndices, noise_constr_dz, sQ(timeStep),
					   sQp(timeStep), sZp(timeStep)));

			factorsByTime.emplace(
				t, boost::make_shared<FactorConstraintsAccIndep>(
					   aMBS, indepCoordIndices, noise_constr_dz, sQ(timeStep),
					   sQp(timeStep), sQpp(timeStep), sZpp(timeStep)));
		}

		// if (timeStep < N - 1)
		{
			factorsByTime.emplace(
				t + dt, boost::make_shared<FactorConstraintsIndep>(
							aMBS, indepCoordIndices, noise_constr_z,
							sZ(timeStep + 1), sQ(timeStep + 1)));

			factorsByTime.emplace(
				t + dt,
				boost::make_shared<FactorConstraintsVelIndep>(
					aMBS, indepCoordIndices, noise_constr_dz, sQ(timeStep + 1),
					sQp(timeStep + 1), sZp(timeStep + 1)));

			factorsByTime.emplace(
				t + dt,
				boost::make_shared<FactorConstraintsAccIndep>(
					aMBS, indepCoordIndices, noise_constr_dz, sQ(timeStep + 1),
					sQp(timeStep + 1), sQpp(timeStep + 1), sZpp(timeStep + 1)));
		}

		// Create initial estimates:
		if (timeStep > 0)
		{
			last_ddq = wholeValues.at<state_t>(sQpp(timeStep - 1));
			last_ddz = wholeValues.at<state_t>(sZpp(timeStep - 1));
		}

		wholeValues.insert(sQpp(timeStep + 1), last_ddq);
		new_timestamps[sQpp(timeStep + 1)] = (timeStep + 1) * dt;

		wholeValues.insert(sZpp(timeStep + 1), last_ddz);
		new_timestamps[sZpp(timeStep + 1)] = (timeStep + 1) * dt;

		// Create initial estimates (so we can run the optimizer)
		wholeValues.insert(sQ(timeStep + 1), last_q);
		new_timestamps[sQ(timeStep + 1)] = (timeStep + 1) * dt;

		wholeValues.insert(sZ(timeStep + 1), last_z);
		new_timestamps[sZ(timeStep + 1)] = (timeStep + 1) * dt;

		wholeValues.insert(sQp(timeStep + 1), last_dq);
		new_timestamps[sQp(timeStep + 1)] = (timeStep + 1) * dt;

		wholeValues.insert(sZp(timeStep + 1), last_dz);
		new_timestamps[sZp(timeStep + 1)] = (timeStep + 1) * dt;

		keysByTime.emplace(t, sQ(timeStep));
		keysByTime.emplace(t, sQp(timeStep));
		keysByTime.emplace(t, sQpp(timeStep));
		keysByTime.emplace(t, sQ(timeStep + 1));
		keysByTime.emplace(t, sQp(timeStep + 1));
		keysByTime.emplace(t, sQpp(timeStep + 1));

		keysByTime.emplace(t, sZ(timeStep));
		keysByTime.emplace(t, sZp(timeStep));
		keysByTime.emplace(t, sZpp(timeStep));
		keysByTime.emplace(t, sZ(timeStep + 1));
		keysByTime.emplace(t, sZp(timeStep + 1));
		keysByTime.emplace(t, sZpp(timeStep + 1));

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

		// For any timestep > 0, add anchoring prior factors to trust the
		// first q,dq values in the window. Not neccesary for t==0 just
		// before we already have prior factors defined (permanent ones, for
		// t=0).
		if (firstTimeIndexInWindow > 0)
		{
			const state_t z_init_win =
				wholeValues.at<state_t>(sZ(firstTimeIndexInWindow));
			const state_t dz_init_win =
				wholeValues.at<state_t>(sZp(firstTimeIndexInWindow));

			// Create Prior factors:
			fgWindow.emplace_shared<gtsam::PriorFactor<state_t>>(
				sZ(firstTimeIndexInWindow), z_init_win, noise_prior_z_0);
			fgWindow.emplace_shared<gtsam::PriorFactor<state_t>>(
				sZp(firstTimeIndexInWindow), dz_init_win, noise_prior_dz_0);
		}

#if 0
    lp.iterationHook = [&fgWindow](size_t iter, double errBef,
                                   double errAfter) {
      const auto N = fgWindow.size();
      std::cout << "LM iter #" << iter << " rmse: " << std::sqrt(errBef / N)
                << " -> " << std::sqrt(errAfter / N) << std::endl;
    };
#endif

#if 0
    fgWindow.print("\n\n***** FG window input*****\n");
    valuesWindow.print("\n\n ******* valuesWindow ****** \n");
#endif

		gtsam::LevenbergMarquardtOptimizer lm(fgWindow, valuesWindow, lp);
		const auto& estimated = lm.optimize();

#if 0
    // std::cout << " === INIT:\n";
    // fgWindow.printErrors(valuesWindow);
    // valuesWindow.print("valuesWindow");
    std::cout << " === FINAL:\n";

    const auto lambdaPrintError = [](const gtsam::Factor *, double err,
                                     size_t) { return err > 0.1; };

    fgWindow.printErrors(estimated, "", gtsam::DefaultKeyFormatter,
                         lambdaPrintError);
    // estimated.print("estimated:");
#endif

		if (!arg_do_not_show_error_progress.isSet())
		{
			const double errorBeforeLM = fgWindow.error(valuesWindow);
			const double errorAfterLM = fgWindow.error(estimated);
			const auto numFactorsLM = fgWindow.size();

			std::cout << "n=" << timeStep << "/" << N << " sliding window LM"
					  << " before RMSE="
					  << std::sqrt(errorBeforeLM / numFactorsLM)
					  << " after RMSE="
					  << std::sqrt(errorAfterLM / numFactorsLM)
					  << " numFactors=" << numFactorsLM
					  << " iters:" << lm.iterations()
					  << " q= " << last_q.transpose() << "\n";
		}

		ASSERT_(lm.iterations() > 0);

		// save/update the last N values (older are "more refined"):
		for (const auto& kv : estimated)
		{
			// Update list of all values (for final batch):
			if (wholeValues.find(kv.key) == wholeValues.end())
				wholeValues.insert(kv.key, kv.value);
			else
				wholeValues.update(kv.key, estimated.at<state_t>(kv.key));
		}

		// Update values in vectors for saving to disk:
		lambda_Values_toQ_DQ_DDQ(estimated);
	}

	// Build FG with all factors:
	if (buildWholeFG)
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
				  << " RMSE before=" << std::sqrt(errorBeforeLM / numFactorsLM)
				  << " numFactors=" << numFactorsLM << "\n"
				  << " RMSE after=" << std::sqrt(errorAfterLM / numFactorsLM)
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

	Qs.saveToTextFile(
		prefix + "q.txt"s, {}, false, "% TIMESTAMP q[0]  ... q[n]\n");
	dotQs.saveToTextFile(
		prefix + "dq.txt"s, {}, false, "% TIMESTAMP dq[0]  ... dq[n]\n");
	ddotQs.saveToTextFile(
		prefix + "ddq.txt"s, {}, false, "% TIMESTAMP ddq[0]  ... ddq[n]\n");
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
