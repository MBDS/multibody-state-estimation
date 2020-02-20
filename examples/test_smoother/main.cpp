// Example of dynamics
// ------------------------------------------------------------
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/slam/PriorFactor.h>
#include <iostream>
#include <fstream>
#include <sparsembs/CAssembledRigidModel.h>
#include <sparsembs/CModelDefinition.h>
#include <sparsembs/FactorDynamics.h>
#include <sparsembs/FactorConstraints.h>
#include <sparsembs/FactorTrapInt.h>
#include <sparsembs/dynamic-simulators.h>
#include <sparsembs/model-examples.h>

void test_smoother()
{
	using gtsam::symbol_shorthand::A;
	using gtsam::symbol_shorthand::Q;
	using gtsam::symbol_shorthand::V;
	using namespace sparsembs;

	// Create the multibody object:
	CModelDefinition model;
	sparsembs::buildFourBarsMBS(model);

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
	const auto n = aMBS->m_q.size();
	const auto m = aMBS->m_Phi_q.getNumRows();

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
	auto noise_constr_q = gtsam::noiseModel::Isotropic::Sigma(m, 0.01);

	const double dt = 0.005;
	const double t_end = 5.0;
	double t = 0;
	unsigned int N = static_cast<unsigned int>(t_end / dt);

	// Create null vector, for use in velocity and accelerations:
	const state_t zeros = gtsam::Vector(gtsam::Vector::Zero(n, 1));

	// Create a feasible Q(0):
	aMBS->m_q.setZero();
	aMBS->m_dotq.setZero();
	aMBS->m_ddotq.setZero();

	CAssembledRigidModel::TComputeDependentParams cdp;  // default params
	CAssembledRigidModel::TComputeDependentResults cdr;
	// Solve the position problem:
	aMBS->m_q[0] = 1;
	aMBS->m_q[1] = 0.1;
	aMBS->m_q[3] = 5;

	aMBS->computeDependentPosVelAcc(indep_coord_indices, true, true, cdp, cdr);
	std::cout << "Position problem final |Phi(q)|=" << cdr.pos_final_phi
			  << "\n";
	ASSERT_BELOW_(cdr.pos_final_phi, 1e-4);

	// Extract m_q from the assembled multibody problem:
	state_t q_0 = gtsam::Vector(aMBS->m_q);
	std::cout << "q0: " << q_0.transpose() << "\n";
	state_t last_q = q_0;

	// Create Prior factors:
	new_factors.emplace_shared<gtsam::PriorFactor<state_t>>(
		Q(0), q_0, noise_prior_q);
	new_factors.emplace_shared<gtsam::PriorFactor<state_t>>(
		V(0), zeros, noise_prior_dq);

	gtsam::ISAM2Params isam2params;
	isam2params.relinearizeThreshold = 0.001;
	isam2params.cacheLinearizedFactors = false;  // Default=true
	isam2params.evaluateNonlinearError = true;  // for debugging mostly

	gtsam::ISAM2 isam2(isam2params);
	gtsam::Values estimated;

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
		// if ((nn % 2) == 0)
		new_factors.emplace_shared<FactorConstraints>(
			aMBS, noise_constr_q, Q(nn));

		// TODO: Add velocity constraints!!

		// Create initial estimates:
		if (!isam2.valueExists(Q(nn)) && !new_values.exists(Q(nn)))
			new_values.insert(Q(nn), last_q);
		if (!isam2.valueExists(V(nn)) && !new_values.exists(V(nn)))
			new_values.insert(V(nn), zeros);
		if (!isam2.valueExists(A(nn)) && !new_values.exists(A(nn)))
			new_values.insert(A(nn), zeros);

		// Create initial estimates (so we can run the optimizer)
		new_values.insert(Q(nn + 1), last_q);
		new_values.insert(V(nn + 1), zeros);
		new_values.insert(A(nn + 1), zeros);

		// Run iSAM every N steps:
		const int SMOOTHER_RUN_PERIOD = 3;
		if ((SMOOTHER_RUN_PERIOD - 1) == (nn % SMOOTHER_RUN_PERIOD))
		{
			// new_factors.print("new_factors:");

			gtsam::ISAM2UpdateParams updateParams;  // defaults

			gtsam::ISAM2Result res =
				isam2.update(new_factors, new_values, updateParams);

			estimated = isam2.calculateEstimate();
			last_q = estimated.at<state_t>(Q(nn));

			std::cout << "Running smoother at t=" << nn << "/" << N
					  << " errorBefore=" << *res.errorBefore
					  << " errorAfter=" << *res.errorAfter << "\n";

			// reset containers:
			new_factors = gtsam::NonlinearFactorGraph();
			new_values.clear();
		}
	}

	// estimated.print("FINAL VALUES:");

	// Save states to files:
	mrpt::math::CMatrixDouble Qs(N, n), dotQs(N, n), ddotQs(N, n);
	for (unsigned int step = 0; step < N; step++)
	{
		const state_t q_val = estimated.at<state_t>(Q(step));
		const state_t dq_val = estimated.at<state_t>(V(step));
		const state_t ddq_val = estimated.at<state_t>(A(step));

		Qs.row(step) = q_val;
		dotQs.row(step) = dq_val;
		ddotQs.row(step) = ddq_val;
	}
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
