/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// These tests exercise the thread-safety story required to evaluate GTSAM
// factors from several threads at once (e.g. GTSAM built with TBB): each
// factor must own a fully independent, private clone of the
// AssembledRigidModel/solver it was constructed from, so concurrent
// evaluateError() calls on different factors never touch shared mutable
// state.

#include <gtest/gtest.h>

#include <atomic>
#include <thread>
#include <vector>

#include <gtsam/inference/Symbol.h>
#include <mbse/AssembledRigidModel.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/factors/FactorConstraints.h>
#include <mbse/factors/FactorDynamics.h>
#include <mbse/model-examples.h>

using namespace mbse;

// ---------------------------------------------------------------------------
// AssembledRigidModel::clone() must produce a fully independent state,
// sharing only the (immutable) topology.
// ---------------------------------------------------------------------------
TEST(ThreadSafety, ModelCloneIsolation)
{
	const ModelDefinition model = mbse::buildFourBarsMBS();
	std::shared_ptr<AssembledRigidModel> original = model.assembleRigidMBS();

	auto clone = original->clone();
	ASSERT_TRUE(clone != nullptr);
	ASSERT_NE(clone.get(), original.get());

	const Eigen::VectorXd q0 = original->q_;

	// Mutating the clone must not be visible in the original, and vice versa:
	clone->q_.setConstant(clone->q_.size(), 123.0);
	EXPECT_TRUE(original->q_.isApprox(q0));
	EXPECT_FALSE(clone->q_.isApprox(q0));

	original->q_.setConstant(original->q_.size(), -7.0);
	EXPECT_FALSE(clone->q_.isApprox(original->q_));

	// Given the SAME state, both must compute identical Jacobians/constraint
	// values (they were built from the same topology):
	original->q_ = q0;
	clone->q_ = q0;
	original->realize_operating_point();
	clone->realize_operating_point();
	original->update_numeric_Phi_and_Jacobians();
	clone->update_numeric_Phi_and_Jacobians();

	EXPECT_TRUE(original->Phi_.isApprox(clone->Phi_, 1e-12));
	EXPECT_TRUE(
		original->Phi_q_.asDense().isApprox(clone->Phi_q_.asDense(), 1e-12));
}

// ---------------------------------------------------------------------------
// CDynamicSimulatorBase::clone() must produce an independent, already
// prepared solver bound to its own independent model clone.
// ---------------------------------------------------------------------------
TEST(ThreadSafety, DynamicSimulatorCloneIndependence)
{
	const ModelDefinition model = mbse::buildFourBarsMBS();
	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_R_matrix_dense original(aMBS);
	original.prepare();

	CDynamicSimulatorBase::Ptr clonedBase = original.clone();
	ASSERT_TRUE(clonedBase != nullptr);

	// The clone must NOT share the same underlying model instance:
	EXPECT_NE(original.get_model().get(), clonedBase->get_model().get());

	const Eigen::VectorXd q0 = original.get_model()->q_;

	// Perturbing and solving on the clone must not touch the original:
	clonedBase->get_model_non_const()->q_[0] += 0.05;
	Eigen::VectorXd ddq;
	clonedBase->solve_ddotq(0.0, ddq);

	EXPECT_TRUE(original.get_model()->q_.isApprox(q0));

	// The clone must be independently usable (already prepared):
	Eigen::VectorXd ddq_orig;
	EXPECT_NO_THROW(original.solve_ddotq(0.0, ddq_orig));
}

// ---------------------------------------------------------------------------
// Each factor must own a private clone: constructing several FactorDynamics
// from the same template solver/model, then evaluating each of them many
// times concurrently from different threads, must be race-free and always
// reproduce the single-threaded reference result.
// ---------------------------------------------------------------------------
TEST(ThreadSafety, ConcurrentFactorDynamicsEvaluation)
{
	const ModelDefinition model = mbse::buildFourBarsMBS();
	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	// A single "template" solver/model, shared only as a source to clone
	// from -- never touched again after the factors below are constructed.
	CDynamicSimulator_R_matrix_dense templateSolver(aMBS);
	templateSolver.prepare();

	const auto n = aMBS->q_.size();
	auto noise = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);

	constexpr int kNumFactors = 8;
	std::vector<FactorDynamics::shared_ptr> factors;
	std::vector<gtsam::Vector> q_values, dq_values, ddq_values;
	std::vector<gtsam::Vector> expected_errors;

	for (int i = 0; i < kNumFactors; i++)
	{
		factors.push_back(std::make_shared<FactorDynamics>(
			&templateSolver, noise, gtsam::Symbol('q', i),
			gtsam::Symbol('v', i), gtsam::Symbol('a', i)));

		gtsam::Vector q = aMBS->q_;
		q[0] += 0.01 * i;  // give each factor a distinct, valid state
		const gtsam::Vector dq = Eigen::VectorXd::Zero(n);
		const gtsam::Vector ddq = Eigen::VectorXd::Zero(n);

		q_values.push_back(q);
		dq_values.push_back(dq);
		ddq_values.push_back(ddq);
	}

	// Single-threaded reference results:
	for (int i = 0; i < kNumFactors; i++)
		expected_errors.push_back(factors[i]->evaluateError(
			q_values[i], dq_values[i], ddq_values[i]));

	// Now hammer every factor concurrently from its own thread, many times,
	// and check the result never deviates from the single-threaded reference
	// (a deviation would mean cross-talk between factors' internal state).
	std::atomic<bool> mismatch{false};
	constexpr int kIters = 40;

	auto worker = [&](int idx)
	{
		for (int it = 0; it < kIters; it++)
		{
			const auto err = factors[idx]->evaluateError(
				q_values[idx], dq_values[idx], ddq_values[idx]);
			if (!err.isApprox(expected_errors[idx], 1e-9)) mismatch = true;
		}
	};

	std::vector<std::thread> threads;
	threads.reserve(kNumFactors);
	for (int i = 0; i < kNumFactors; i++) threads.emplace_back(worker, i);
	for (auto& t : threads) t.join();

	EXPECT_FALSE(mismatch.load());
}

// ---------------------------------------------------------------------------
// Same idea, but for a kinematic-only factor (FactorConstraints), which owns
// a private AssembledRigidModel clone directly (no solver involved).
// ---------------------------------------------------------------------------
TEST(ThreadSafety, ConcurrentFactorConstraintsEvaluation)
{
	const ModelDefinition model = mbse::buildFourBarsMBS();
	std::shared_ptr<AssembledRigidModel> aMBS = model.assembleRigidMBS();

	const auto m = aMBS->Phi_.size();
	auto noise = gtsam::noiseModel::Isotropic::Sigma(m, 0.1);

	constexpr int kNumFactors = 8;
	std::vector<FactorConstraints::shared_ptr> factors;
	std::vector<gtsam::Vector> q_values;
	std::vector<gtsam::Vector> expected_errors;

	for (int i = 0; i < kNumFactors; i++)
	{
		factors.push_back(std::make_shared<FactorConstraints>(
			aMBS, noise, gtsam::Symbol('q', i)));

		gtsam::Vector q = aMBS->q_;
		q[0] += 0.01 * i;
		q_values.push_back(q);
	}

	for (int i = 0; i < kNumFactors; i++)
		expected_errors.push_back(factors[i]->evaluateError(q_values[i]));

	std::atomic<bool> mismatch{false};
	constexpr int kIters = 40;

	auto worker = [&](int idx)
	{
		for (int it = 0; it < kIters; it++)
		{
			const auto err = factors[idx]->evaluateError(q_values[idx]);
			if (!err.isApprox(expected_errors[idx], 1e-9)) mismatch = true;
		}
	};

	std::vector<std::thread> threads;
	threads.reserve(kNumFactors);
	for (int i = 0; i < kNumFactors; i++) threads.emplace_back(worker, i);
	for (auto& t : threads) t.join();

	EXPECT_FALSE(mismatch.load());

	// The original model, only ever used as a clone template, must remain
	// untouched by any of the concurrent evaluations above:
	EXPECT_TRUE(aMBS->q_.isApprox(model.assembleRigidMBS()->q_));
}
