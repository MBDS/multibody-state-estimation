/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

// Example of numerical integration using factor graphs:
// ------------------------------------------------------------
#include <gtest/gtest.h>

#include <mbse/factors/FactorEulerInt.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/factorTesting.h>
#include <mrpt/core/exceptions.h>  // ASSERT_()
#include <iostream>

using namespace std;

TEST(FactorEulerInt, optimize)
{
	using gtsam::symbol_shorthand::V;
	using gtsam::symbol_shorthand::X;
	using namespace mbse;

	// Create the factor graph:
	gtsam::NonlinearFactorGraph graph;
	gtsam::Values initValues;

	// Add factors:
	auto noise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector2(1, 1));
	auto noise_prior =
		gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector2(1, 1));

	const double dt = 1.0;

	graph.emplace_shared<FactorEulerInt>(dt, noise, X(1), X(2), V(1));
	graph.emplace_shared<FactorEulerInt>(dt, noise, X(2), X(3), V(2));

	const state_t zeros = gtsam::Vector(gtsam::Vector2(0.0, 0.0));

	const state_t prior_x1 = gtsam::Vector(gtsam::Vector2(0.0, 0.0));
	const state_t prior_v1 = gtsam::Vector(gtsam::Vector2(1.0, 2.0));
	const state_t prior_v2 = gtsam::Vector(gtsam::Vector2(2.0, 3.0));

	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		X(1), prior_x1, noise_prior);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(1), prior_v1, noise_prior);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(2), prior_v2, noise_prior);

	initValues.insert(X(1), zeros);
	initValues.insert(X(2), zeros);
	initValues.insert(X(3), zeros);
	initValues.insert(V(1), zeros);
	initValues.insert(V(2), zeros);

	// Run optimizer:
	// graph.print("Factor graph: ");
	// initValues.print("initValues: ");
	gtsam::LevenbergMarquardtOptimizer optimizer(graph, initValues);

	const auto& optimValues = optimizer.optimize();

	// Process results:
	// optimValues.print("optimValues");
	/* Expected values:
	 * optimValuesValues with 5 values:
	 * Value v1: (N9sparsembs7state_tE) 1 2
	 * Value v2: (N9sparsembs7state_tE) 2 3
	 * Value x1: (N9sparsembs7state_tE) -3.79966e-10 -6.59941e-10
	 * Value x2: (N9sparsembs7state_tE) 1 2
	 * Value x3: (N9sparsembs7state_tE) 3 5
	 */

	EXPECT_NEAR(optimValues.at<state_t>(V(1))[0], 1.0, 1e-6);
	EXPECT_NEAR(optimValues.at<state_t>(V(1))[1], 2.0, 1e-6);

	EXPECT_NEAR(optimValues.at<state_t>(V(2))[0], 2.0, 1e-6);
	EXPECT_NEAR(optimValues.at<state_t>(V(2))[1], 3.0, 1e-6);

	EXPECT_NEAR(optimValues.at<state_t>(X(1))[0], 0.0, 1e-6);
	EXPECT_NEAR(optimValues.at<state_t>(X(1))[1], 0.0, 1e-6);

	EXPECT_NEAR(optimValues.at<state_t>(X(2))[0], 1.0, 1e-6);
	EXPECT_NEAR(optimValues.at<state_t>(X(2))[1], 2.0, 1e-6);

	EXPECT_NEAR(optimValues.at<state_t>(X(3))[0], 3.0, 1e-6);
	EXPECT_NEAR(optimValues.at<state_t>(X(3))[1], 5.0, 1e-6);
}

TEST(FactorEulerInt, Jacobian)
{
// Use MRPT macros for testing as "EXPECT"
#define EXPECT ASSERT_
	const std::string name_ = "FactorEulerInt";

	using gtsam::symbol_shorthand::V;
	using gtsam::symbol_shorthand::X;
	using namespace mbse;

	auto noise = gtsam::noiseModel::Diagonal::Sigmas(gtsam::Vector2(1, 1));
	const double dt = 1e-3;

	FactorEulerInt factor(dt, noise, X(1), X(2), V(1));

	// Set the linearization point
	gtsam::Values values;

	const state_t x1 = gtsam::Vector(gtsam::Vector2(1.0, 2.0));
	const state_t x2 = gtsam::Vector(gtsam::Vector2(3.0, 4.0));
	const state_t v1 = gtsam::Vector(gtsam::Vector2(5.0, 6.0));

	values.insert(X(1), x1);
	values.insert(X(2), x2);
	values.insert(V(1), v1);

	EXPECT_CORRECT_FACTOR_JACOBIANS(
		factor, values, 1e-7 /*diff*/, 1e-6 /*tolerance*/);
}
