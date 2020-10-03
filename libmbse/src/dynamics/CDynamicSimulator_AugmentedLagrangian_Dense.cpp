/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CAssembledRigidModel.h>
#include <mbse/dynamics/dynamic-simulators.h>

using namespace mbse;
using namespace Eigen;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense solver with the Augmented Lagrange formulation (ALF)
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_AugmentedLagrangian_Dense::
	CDynamicSimulator_AugmentedLagrangian_Dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBasePenalty(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_AugmentedLagrangian_Dense::internal_prepare()
{
	timelog().enter("solver_prepare");

	arm_->buildMassMatrix_dense(M_);
	M_ldlt_.compute(M_);

	timelog().leave("solver_prepare");
}

CDynamicSimulator_AugmentedLagrangian_Dense::
	~CDynamicSimulator_AugmentedLagrangian_Dense()
{
}

void CDynamicSimulator_AugmentedLagrangian_Dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	const size_t nDepCoords = arm_->q_.size();

	if (lagrangre)
		throw std::runtime_error(
			"This class can't solve for lagrange multipliers!");

	timelog().enter("solver_ddotq");

	// Iterative solution to the Augmented Lagrangian Formulation (ALF):
	// ---------------------------------------------------------------------

	// 1) M \ddot{q}_0 = Q
	// ---------------------------
	// Get "Q":
	// KLU leaves solution in the same place than the input RHS vector:
	Eigen::VectorXd ddotq_prev(nDepCoords), ddotq_next(nDepCoords);
	this->build_RHS(&ddotq_prev[0] /* Q */, NULL /* we don't need "c" */);

	ddotq_prev = M_ldlt_.solve(ddotq_prev);

	// 2) Iterate:
	// ---------------------------
	//
	// [ M + alpha * Phi_q^t * Phi_q ] \ddot{q}_i+1 = RHS
	//
	// \--------------v-------------/
	//                =A
	//
	// RHS = M*\ddot{q}_i - alpha * Phi_q^t* [ \dot{Phi}_q * \dot{q} + 2 * xi *
	// omega * \dot{q} + omega^2 * Phi  ]
	//

	// Update numeric values of the constraint Jacobians:
	arm_->update_numeric_Phi_and_Jacobians();

	arm_->getPhi_q_dense(Phi_q_);
	A_ = M_ + params_penalty.alpha * Phi_q_.transpose() * Phi_q_;

	A_lu_.compute(A_);

	// Build the RHS vector:
	// RHS = M*\ddot{q}_i -  Phi_q^t* alpha * [ \dot{Phi}_q * \dot{q} + 2 * xi *
	// omega * \dot{Phi} + omega^2 * Phi  ]
	//                               \ ------------------------------------v
	//                               --------------------------------------/
	//                                                                    = b
	timelog().enter("solver_ddotq.build_rhs");

	arm_->dotPhi_q_.asDense(dotPhi_q_);

	const Eigen::MatrixXd RHS2 =
		params_penalty.alpha * Phi_q_.transpose() *
		(dotPhi_q_ * arm_->dotq_ +
		 2 * params_penalty.xi * params_penalty.w * arm_->dotPhi_ +
		 params_penalty.w * params_penalty.w * arm_->Phi_);

	timelog().leave("solver_ddotq.build_rhs");

	// Solve linear system:
	// -----------------------------------
	timelog().enter("solver_ddotq.solve");

	Eigen::VectorXd RHS(nDepCoords);

	const double MAX_DDOT_INCR_NORM = 1e-4 * nDepCoords;
	const size_t MAX_ITERS = 10;

	double ddot_incr_norm;
	size_t iter = 0;
	do
	{
		// RHS = M*\ddot{q}_i - RHS2
		RHS = (M_ * ddotq_prev) - RHS2;
		ddotq_next = A_lu_.solve(RHS);

		ddot_incr_norm = (ddotq_next - ddotq_prev).norm();
		// cout << "iter: " << iter<< endl << "prev: " << ddotq_prev.transpose()
		// << "\nnext: " << ddotq_next.transpose() << "\n  norm: " <<
		// ddot_incr_norm << endl << endl;

		ddotq_prev = ddotq_next;
	} while (ddot_incr_norm > MAX_DDOT_INCR_NORM && ++iter < MAX_ITERS);

	ddot_q.swap(ddotq_next);

	timelog().leave("solver_ddotq.solve");

	ASSERTDEBMSG_(
		((RHS.array() == RHS.array()).all()), "NaN found in result ddotq");

	timelog().leave("solver_ddotq");
}

/** Integrators will call this after each time step */
void CDynamicSimulator_AugmentedLagrangian_Dense::post_iteration(double t)
{
#if 0

	const size_t nConstraints = arm_->Phi_.size();

	// ---------------------------
	// Projection in q
	// ---------------------------
	const Eigen::VectorXd q0 = arm_->q_;
	Eigen::VectorXd Lambda(nConstraints);
	Lambda.setZero();

	Eigen::VectorXd rhs, Aq, Aqp;

	for (int i=0;i<3;i++)
	{
		// Update numeric values of the constraint Jacobians:
		arm_->update_numeric_Phi_and_Jacobians();
		arm_->getPhi_q_dense(Phi_q_);

		Lambda += params_penalty.alpha * arm_->Phi_;

		// Solve for increments Aq:
		//
		// -> Lambda = Lambda + alpha * Phi
		// -> [M+alpha * Phi_q^t * Phi_q] Aq = -[ M (qi-q0) + Phi_q^t * Lambda ]
		rhs = M_ *( q0 - arm_->q_ ) - Phi_q_.transpose() * Lambda;

		Aq = A_lu_.solve(rhs);
		arm_->q_ += Aq;

		cout << "iter: " << i << " |Aq|=" << Aq.norm() << endl;
	}
	cout << "\n";

	// ---------------------------
	// Projection in dot{q}
	// ---------------------------
//	Eigen::FullPivLU<Eigen::MatrixXd> lu;
//	lu.compute(Phi_q_);
//	const Eigen::MatrixXd V = lu.kernel();
//	arm_->dotq_ = (V*V.transpose()) * arm_->dotq_;

	timelog().leave("solver.post_iteration");
#endif  // 0	timelog().enter("solver.post_iteration");
}
