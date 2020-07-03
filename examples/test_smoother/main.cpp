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
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam_unstable/nonlinear/BatchFixedLagSmoother.h>
#include <gtsam/slam/PriorFactor.h>
#include <iostream>
#include <fstream>
#include <mbse/CAssembledRigidModel.h>
#include <mbse/CModelDefinition.h>
#include <mbse/factors/FactorDynamics.h>
#include <mbse/factors/FactorConstraints.h>
#include <mbse/factors/FactorConstraintsVel.h>
#include <mbse/factors/FactorEulerInt.h>
#include <mbse/factors/FactorTrapInt.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/model-examples.h>

void test_smoother()
{
	using gtsam::symbol_shorthand::A;
	using gtsam::symbol_shorthand::Q;
	using gtsam::symbol_shorthand::V;
	using namespace mbse;

	// Create the multibody object:
	CModelDefinition model;
	mbse::buildFourBarsMBS(model);

	std::shared_ptr<CAssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_R_matrix_dense dynSimul(aMBS);

	// Must be called before solve_ddotq(), needed inside the dynamics factors
	dynSimul.prepare();

	// iSAM2: Incrementally build the graph:
	gtsam::NonlinearFactorGraph new_factors;
	gtsam::Values new_values;

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->q_.size();
	const auto m = aMBS->Phi_q_.getNumRows();

	const double noise_vel_sigma = 0.01, noise_acc_sigma = 0.01;

	auto noise_vel = gtsam::noiseModel::Isotropic::Sigma(n, noise_vel_sigma);
	auto noise_acc = gtsam::noiseModel::Isotropic::Sigma(n, noise_acc_sigma);

	// x1, *y1*, x2, y2
	// 0   1     2   3
	std::vector<size_t> indep_coord_indices;
	indep_coord_indices.push_back(0);

	// Velocity prior: large sigma for all dq(i), except dq(i_indep)
	gtsam::Vector prior_dq_sigmas;
	const double large_std = 1e6;
	const double small_std = 1e-3;
	prior_dq_sigmas.setConstant(n, large_std);
	for (auto idx : indep_coord_indices) prior_dq_sigmas(idx) = small_std;

	auto noise_prior_dq = gtsam::noiseModel::Diagonal::Sigmas(prior_dq_sigmas);
	auto noise_prior_q = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);
	auto noise_dyn = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);
	auto noise_constr_q = gtsam::noiseModel::Isotropic::Sigma(m, 0.1);
	auto noise_constr_dq = gtsam::noiseModel::Isotropic::Sigma(m, 0.1);

	const double dt = 0.005;
	const double t_end = 5.0;
	double t = 0;
	unsigned int N = static_cast<unsigned int>(t_end / dt);

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Create a feasible Q(0):
	aMBS->q_.setZero();
	aMBS->dotq_.setZero();
	aMBS->ddotq_.setZero();

	CAssembledRigidModel::TComputeDependentParams cdp;  // default params
	CAssembledRigidModel::TComputeDependentResults cdr;
	// Solve the position problem:
	aMBS->q_[0] = 1;
	aMBS->q_[1] = 0.1;
	aMBS->q_[3] = 5;

	aMBS->computeDependentPosVelAcc(indep_coord_indices, true, true, cdp, cdr);
	std::cout << "Position problem final |Phi(q)|=" << cdr.pos_final_phi
			  << "\n";
	ASSERT_BELOW_(cdr.pos_final_phi, 1e-4);

	// Extract q_ from the assembled multibody problem:
	state_t q_0 = gtsam::Vector(aMBS->q_);
	std::cout << "q0: " << q_0.transpose() << "\n";
	state_t last_q = q_0, last_dq = zeros, last_ddq = zeros;

	// Create Prior factors:
	new_factors.emplace_shared<gtsam::PriorFactor<state_t>>(
		Q(0), q_0, noise_prior_q);
	new_factors.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(0), zeros, noise_prior_dq);

	const double lag = 0.1;  // seconds
	gtsam::BatchFixedLagSmoother smoother(lag);

	smoother.params().maxIterations = 15;

	gtsam::Values estimated;
	gtsam::FixedLagSmoother::KeyTimestampMap new_timestamps;

	new_values.insert(Q(0), last_q);
	new_values.insert(V(0), last_dq);
	new_values.insert(A(0), last_ddq);
	new_timestamps[Q(0)] = 0 * dt;
	new_timestamps[V(0)] = 0 * dt;
	new_timestamps[A(0)] = 0 * dt;

	// Save states to files:
	mrpt::math::CMatrixDouble Qs(N + 1, n), dotQs(N + 1, n), ddotQs(N + 1, n);

	for (unsigned int nn = 0; nn < N; nn++, t += dt)
	{
		// Create Trapezoidal Integrator factors:
		new_factors.emplace_shared<FactorTrapInt>(
			dt, noise_vel, Q(nn), Q(nn + 1), V(nn), V(nn + 1));
		new_factors.emplace_shared<FactorTrapInt>(
			dt, noise_acc, V(nn), V(nn + 1), A(nn), A(nn + 1));

		// Create Dynamics factors:
		new_factors.emplace_shared<FactorDynamics>(
			&dynSimul, noise_dyn, Q(nn), V(nn), A(nn));

		// Add dependent-coordinates constraint factor:
		new_factors.emplace_shared<FactorConstraints>(
			aMBS, noise_constr_q, Q(nn));
		new_factors.emplace_shared<FactorConstraintsVel>(
			aMBS, noise_constr_dq, Q(nn), V(nn));

		// Create initial estimates:
		new_values.insert(A(nn + 1), last_ddq);
		new_timestamps[A(nn + 1)] = (nn + 1) * dt;

		// Create initial estimates (so we can run the optimizer)
		new_values.insert(Q(nn + 1), last_q);
		new_timestamps[Q(nn + 1)] = (nn + 1) * dt;

		new_values.insert(V(nn + 1), last_dq);
		new_timestamps[V(nn + 1)] = (nn + 1) * dt;

		// new_values.insert(A(nn + 1), zeros); // not for Euler integrator

		// Run iSAM every N steps:
		// if ((nn % 10) == 0 || nn == (N - 1))  <------------- Commented
		{
			std::cout << "n=" << nn << "/" << N
					  << " keys: " << smoother.timestamps().size() << "\n";
			const auto res =
				smoother.update(new_factors, new_values, new_timestamps);

			std::cout << "Error: " << res.getError() << "\n";

			// reset containers:
			new_factors = gtsam::NonlinearFactorGraph();
			new_values.clear();
			new_timestamps.clear();

			estimated = smoother.calculateEstimate();

			// save/update the last N values (older are "more refined"):
			for (const auto& kv : estimated)
			{
				gtsam::Symbol s(kv.key);
				const auto step = s.index();
				const state_t val = estimated.at<state_t>(s);
				switch (s.chr())
				{
					case 'q':
						Qs.row(step) = val;
						break;
					case 'v':
						dotQs.row(step) = val;
						break;
					case 'a':
						ddotQs.row(step) = val;
						break;
				};
			}
		}
	}

	// new_factors.print("new_factors:");
	//	const double errorBefore = smoother.getFactors().error(estimated);

	const double errorAfter = smoother.getFactors().error(estimated);

	//	last_q = estimated.at<state_t>(Q(nn));
	//	last_dq = estimated.at<state_t>(V(nn));
	//	last_ddq = estimated.at<state_t>(A(nn));

	const auto numFactors = smoother.getFactors().size();

	std::cout << "ErrorAfter=" << errorAfter
			  << " RMSE=" << std::sqrt(errorAfter / numFactors)
			  << " numFactors=" << numFactors << "\n";

	// estimated.print("FINAL VALUES:");

	/*for (const auto& key_timestamp : smoother.timestamps())
	{
		// std::cout << "Key: " << key_timestamp.first
		// << " Time: " << key_timestamp.second << std::endl;

		const gtsam::Symbol s(key_timestamp.first);


	}*/

	std::cout << "Saving results to TXT files...\n";
	Qs.saveToTextFile("q.txt");
	dotQs.saveToTextFile("dq.txt");
	ddotQs.saveToTextFile("ddq.txt");
}

int main()
{
	try
	{
		test_smoother();
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
	}
}
