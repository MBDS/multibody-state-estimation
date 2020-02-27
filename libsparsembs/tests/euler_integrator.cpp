// Example of numerical integration using factor graphs:
// ------------------------------------------------------------
#include <gtest/gtest.h>

#include <sparsembs/FactorEulerInt.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <iostream>

using namespace std;

TEST(numerical_integrators, euler)
{
	using gtsam::symbol_shorthand::V;
	using gtsam::symbol_shorthand::X;
	using namespace sparsembs;

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
	graph.print("Factor graph: ");
	initValues.print("initValues: ");
	gtsam::LevenbergMarquardtOptimizer optimizer(graph, initValues);

	const auto& optimValues = optimizer.optimize();

	// Process results:
	optimValues.print("optimValues");
}
