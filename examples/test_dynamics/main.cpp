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
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = -0.05;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = std::sqrt(2.0 * 2.0 + 3.0 * 3.0);
		b.mass() = 4;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
}

void test_dynamics()
{
	using gtsam::symbol_shorthand::A;
	using gtsam::symbol_shorthand::Q;
	using gtsam::symbol_shorthand::V;
	// using namespace sparsembs; /Already used in the global section above

	// Create the multibody object:
	CModelDefinition model;
	buildFourBarsMBS(model);

	std::shared_ptr<CAssembledRigidModel> aMBS = model.assembleRigidMBS();
	aMBS->setGravityVector(0, -9.81, 0);

	CDynamicSimulator_R_matrix_dense dynSimul(aMBS);

	// Must be called before solve_ddotq(), needed inside the dynamics factors
	dynSimul.prepare();

	// Create the empty factor graph:
	gtsam::NonlinearFactorGraph graph;
	gtsam::Values initValues;

	// Add factors:
	// Create factor noises:
	const auto n = aMBS->m_q.size();

	const double noise_vel_sigma = 0.01, noise_acc_sigma = 0.01;

	auto noise_vel = gtsam::noiseModel::Isotropic::Sigma(n, noise_vel_sigma);
	auto noise_acc = gtsam::noiseModel::Isotropic::Sigma(n, noise_acc_sigma);

	// x1, *y1*, x2, y2
	// 0   1     2   3
	std::vector<size_t> indep_coord_indices;
	indep_coord_indices.push_back(1);

	// Velocity prior: large sigma for all dq(i), except dq(i_indep)
	gtsam::Vector prior_dq_sigmas;
	const double large_std = 1e6;
	const double small_std = 1e-3;
	prior_dq_sigmas.setConstant(n, large_std);
	for (auto idx : indep_coord_indices) prior_dq_sigmas(idx) = small_std;

	auto noise_prior_dq = gtsam::noiseModel::Diagonal::Sigmas(prior_dq_sigmas);
	auto noise_prior_q = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);
	auto noise_dyn = gtsam::noiseModel::Isotropic::Sigma(n, 0.1);

	const double dt = 0.01;

	// Create Trapezoidal Integrator factors:
	graph.emplace_shared<FactorTrapInt>(dt, noise_vel, Q(0), Q(1), V(0), V(1));
	graph.emplace_shared<FactorTrapInt>(dt, noise_acc, V(0), V(1), A(0), A(1));

	graph.emplace_shared<FactorTrapInt>(dt, noise_vel, Q(1), Q(2), V(1), V(2));
	graph.emplace_shared<FactorTrapInt>(dt, noise_acc, V(1), V(2), A(1), A(2));

	// Create Dynamics factors:
	graph.emplace_shared<FactorDynamics>(
		&dynSimul, noise_dyn, Q(0), V(0), A(0));
	graph.emplace_shared<FactorDynamics>(
		&dynSimul, noise_dyn, Q(1), V(1), A(1));
	graph.emplace_shared<FactorDynamics>(
		&dynSimul, noise_dyn, Q(2), V(2), A(2));

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Create a feasible Q(0):
	aMBS->m_q.setZero();
	aMBS->m_dotq.setZero();
	aMBS->m_ddotq.setZero();

	CAssembledRigidModel::TComputeDependentParams cdp;  // default params
	CAssembledRigidModel::TComputeDependentResults cdr;
	// Solve the position problem:
	aMBS->m_q[1] = 3;
	aMBS->computeDependentPosVelAcc(indep_coord_indices, true, true, cdp, cdr);
	std::cout << "Position problem final |Phi(q)|=" << cdr.pos_final_phi
			  << "\n";

	// Extract m_q from the assembled multibody problem:
	state_t q_0 = gtsam::Vector(aMBS->m_q);
	std::cout << "q0: " << q_0.transpose() << "\n";

	// Create Prior factors:
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(Q(0), q_0, noise_prior_q);
	graph.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(0), zeros, noise_prior_dq);

	// Create initial estimates:
	initValues.insert(Q(0), q_0);
	initValues.insert(Q(1), q_0);
	initValues.insert(Q(2), q_0);
	initValues.insert(V(0), zeros);
	initValues.insert(V(1), zeros);
	initValues.insert(V(2), zeros);
	initValues.insert(A(0), zeros);
	initValues.insert(A(1), zeros);
	initValues.insert(A(2), zeros);

	// Run optimizer:
	std::cout.precision(3);

	graph.print("Factor graph: ");
	initValues.print("initValues: ");
	gtsam::LevenbergMarquardtOptimizer optimizer(graph, initValues);

	const auto& optimValues = optimizer.optimize();

	// Process results:
	optimValues.print("optimValues");

	std::cout << "Initial factors error: " << graph.error(initValues) << "\n";
	std::cout << "Final factors error: " << graph.error(optimValues) << "\n";
	std::cout << "Optimization iterations: " << optimizer.iterations() << "\n";
}

int main()
{
	try
	{
		test_dynamics();
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
	}
}
