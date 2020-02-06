// Example of dynamics
// ------------------------------------------------------------
#include <sparsembs/FactorTrapInt.h>
#include <sparsembs/FactorDynamics.h>
#include <sparsembs/CModelDefinition.h>
#include <sparsembs/CAssembledModelRigid.h>
#include <sparsembs/dynamic-simulators.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <iostream>

using namespace std;
using namespace sparsembs;

void buildFourBarsMBS(CModelDefinition& model)
{
	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 0));
	model.setPointCoords(2, TPoint2D(1, 2));
	model.setPointCoords(3, TPoint2D(4, 0), true /*is fixed*/);

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 1;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = -0.05;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = std::sqrt(2.0 * 2.0 + 3.0 * 3.0);
		b.mass() = 4;
		b.I0() = (1. / 3.) * b.mass() * square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
}

void test_dynamics()
{
	using gtsam::symbol_shorthand::Q;
	using gtsam::symbol_shorthand::V;
	using gtsam::symbol_shorthand::A;
	using namespace sparsembs;

	// Create the multibody object:
	CModelDefinition model;
	buildFourBarsMBS(model);
	
	std::shared_ptr<CAssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_R_matrix_dense dynSimul(aMBS);

	// Create the factor graph:
	gtsam::NonlinearFactorGraph graph;
	gtsam::Values initValues;

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->m_q.size();

	const double noise_vel_sigma = 0.1, noise_acc_sigma = 0.1;

	auto noise_vel = gtsam::noiseModel::Isotropic(n, noise_vel_sigma);
	auto noise_acc = gtsam::noiseModel::Isotropic(n, noise_acc_sigma);

	auto noise_prior_q = gtsam::noiseModel::Isotropic(n, 0.1);
	auto noise_prior_dq = gtsam::noiseModel::Isotropic(n, 0.1);

	auto noise_dyn = gtsam::noiseModel::Isotropic(n, 0.1);

	const double dt = 0.1;

	// Create Trapezoidal Integrator factors:
	graph.emplace_shared<FactorTrapInt>(dt, noise_vel, Q(0), Q(1), V(0), V(1));
	graph.emplace_shared<FactorTrapInt>(dt, noise_acc, V(0), V(1), A(0), A(1));

	graph.emplace_shared<FactorTrapInt>(dt, noise_vel, Q(1), Q(2), V(1), V(2));
	graph.emplace_shared<FactorTrapInt>(dt, noise_acc, V(1), V(2), A(1), A(2));

	// Dynamics:
	graph.emplace_shared<FactorDynamic>(&dynSimul, noise_dyn, Q(0), V(0), A(0));
	graph.emplace_shared<FactorDynamic>(&dynSimul, noise_dyn, Q(1), V(1), A(1));

	// Continue here!

	// Create null vector:
	const state_t zeros = gtsam::Vector(gtsam::Vector2(0.0, 0.0));

	// Create Prior vectors:
	const state_t prior_x1 = gtsam::Vector(gtsam::Vector2(0.0, 0.0));
	const state_t prior_v1 = gtsam::Vector(gtsam::Vector2(1.0, 2.0));
	const state_t prior_v2 = gtsam::Vector(gtsam::Vector2(2.0, 3.0));
	const state_t prior_v3 = gtsam::Vector(gtsam::Vector2(2.0, 3.0));

	// Create Prior factors:
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		X(1), prior_x1, noise_prior);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(1), prior_v1, noise_prior);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(2), prior_v2, noise_prior);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(3), prior_v3, noise_prior);


	// Create initial estimates:
	initValues.insert(X(1), zeros);
	initValues.insert(X(2), zeros);
	initValues.insert(X(3), zeros);
	initValues.insert(V(1), zeros);
	initValues.insert(V(2), zeros);
	initValues.insert(V(3), zeros);

	// Run optimizer:
	graph.print("Factor graph: ");
	initValues.print("initValues: ");
	gtsam::LevenbergMarquardtOptimizer optimizer(graph, initValues);

	const auto& optimValues = optimizer.optimize();

	// Process results:
	optimValues.print("optimValues");
}

int main()
{
	try
	{
		test_dynamics();
	}
	catch (const std::exception & e)
	{
		std::cerr << "Error: " << e.what() << "\n";
	}
}
