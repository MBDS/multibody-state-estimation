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
//  Solver: Dense solver with the index-3 Augmented Lagrange formulation (ALF)
//  with projection of velocities and acelerations
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_ALi3_Dense::CDynamicSimulator_ALi3_Dense(
	const CAssembledRigidModel::Ptr arm_ptr)
	: CDynamicSimulatorBasePenalty(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_ALi3_Dense::internal_prepare()
{
	timelog.enter("solver_prepare");

	arm_->buildMassMatrix_dense(M_);
	M_ldlt_.compute(M_);

	const size_t nConstraints = arm_->Phi_.size();
	Lambda_.setZero(nConstraints);

	timelog.leave("solver_prepare");
}

CDynamicSimulator_ALi3_Dense::~CDynamicSimulator_ALi3_Dense() {}

/** Implement a especific combination of dynamic formulation + integrator.
 *  \return false if it's not implemented, so it should fallback to generic
 * integrator + internal_solve_ddotq()
 */
bool CDynamicSimulator_ALi3_Dense::internal_integrate(
	double t, double dt, const ODE_integrator_t integr)
{
	if (integr != ODE_Trapezoidal) return false;

	const size_t nDepCoords = arm_->q_.size();

	timelog.enter("internal_integrate");

	Eigen::VectorXd Q(nDepCoords);

	const double dt2 = dt * dt;

	const Eigen::VectorXd qp_g = -(2. / dt * arm_->q_ + arm_->dotq_);
	const Eigen::VectorXd qpp_g =
		-(4. / dt2 * arm_->q_ + 4. / dt * arm_->dotq_ + arm_->ddotq_);

	arm_->q_ += dt * arm_->dotq_ + 0.5 * dt * dt * arm_->ddotq_;

	arm_->dotq_ = (2. / dt) * arm_->q_ + qp_g;
	arm_->ddotq_ = (4. / dt2) * arm_->q_ + qpp_g;

	double err = 1;
	int iter = 0;

	const double tol_dyn = 1e-6;
	const int iter_max = 20;

	arm_->update_numeric_Phi_and_Jacobians();
	arm_->getPhi_q_dense(Phi_q_);

	while (err > tol_dyn && iter < iter_max)
	{
		iter++;

		// phi_0 = phi(q,l,x);
		// Get "Q" (may be dynamic)
		this->build_RHS(&Q[0] /* Q */, NULL /* we don't need "c" */);
		// Q = Qg +[0;0;0;0;-k_m*(q(5)-L_0)-c_m*qp(5)];

		Eigen::VectorXd RHS =
			0.25 * dt2 *
			(M_ * arm_->ddotq_ +
			 Phi_q_.transpose() * params_penalty.alpha * arm_->Phi_ +
			 Phi_q_.transpose() * Lambda_ - Q);
		//[K,C]=evalKC(k_m, c_m);

		// f_q = M + 0.5*dt*C+0.25*dt^2*(jac'*alpha*jac+K);
		A_ = M_ +
			  0.25 * dt2 * params_penalty.alpha * Phi_q_.transpose() * Phi_q_;
		A_lu_.compute(A_);

		const Eigen::VectorXd Aq = -A_lu_.solve(RHS);

		arm_->q_ += Aq;
		arm_->dotq_ = (2. / dt) * arm_->q_ + qp_g;
		arm_->ddotq_ = (4. / dt2) * arm_->q_ + qpp_g;

		// phi_0 = phi(q,l,x);
		arm_->update_numeric_Phi_and_Jacobians();
		arm_->getPhi_q_dense(Phi_q_);

		Lambda_ += params_penalty.alpha * arm_->Phi_;
		err = Aq.norm();
	}

	// cout << "iter: " << iter << endl;

	// Proyecciones en velocidad y aceleración (faltan los términos dependientes
	// del timepo, porque en este problema no hay restricciones que dependan
	// explícitamente del tiempo).
	// qp_out = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qp);
	arm_->dotq_ = A_lu_.solve(M_ * arm_->dotq_);

	// phiqpqp_0 = phiqpqp(q, qp, l);
	arm_->dotPhi_q_.asDense(dotPhi_q_);

	// qpp_out = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qpp -
	// 0.25*dt^2*jac'*alpha*phiqpqp_0);
	arm_->ddotq_ = A_lu_.solve(
		M_ * arm_->ddotq_ - 0.25 * dt2 * params_penalty.alpha *
								   Phi_q_.transpose() * dotPhi_q_ *
								   arm_->dotq_);

	timelog.leave("internal_integrate");

	return true;
}

void CDynamicSimulator_ALi3_Dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	const size_t nDepCoords = arm_->q_.size();

	if (lagrangre)
		throw std::runtime_error(
			"This class can't solve for lagrange multipliers!");

	timelog.enter("solver_ddotq");

	// Iterative solution to the Augmented Lagrangian Formulation (ALF):
	// ---------------------------------------------------------------------

	// 1) M \ddot{q}_0 = Q
	// ---------------------------
	// Get "Q":
	// KLU leaves solution in the same place than the input RHS vector:
	Eigen::VectorXd Q(nDepCoords);
	this->build_RHS(&Q[0] /* Q */, NULL /* we don't need "c" */);

	// 2) Iterate:
	// ---------------------------
	//
	// [ M + alpha * Phi_q^t * Phi_q ] \ddot{q}_i+1 = RHS
	//
	// \--------------v-------------/
	//                =A
	//
	// RHS = Q(q,dq) - alpha * Phi_q^t* [ \dot{Phi}_q * \dot{q} + 2 * xi * omega
	// * \dot{Phi} + omega^2 * Phi ] - Phi_q^t * \lambda
	//

	// Update numeric values of the constraint Jacobians:
	arm_->update_numeric_Phi_and_Jacobians();

	arm_->getPhi_q_dense(Phi_q_);
	A_ = M_ + params_penalty.alpha * Phi_q_.transpose() * Phi_q_;

	A_lu_.compute(A_);

	// Build the RHS vector:
	// RHS = Q(q,dq) - alpha * Phi_q^t* [ \dot{Phi}_q * \dot{q} + 2 * xi * omega
	// * \dot{Phi} + omega^2 * Phi ] - Phi_q^t * \lambda
	//

	arm_->dotPhi_q_.asDense(dotPhi_q_);

	// Solve linear system:
	// -----------------------------------
	timelog.enter("solver_ddotq.solve");

	Eigen::VectorXd ddotq_next, ddotq_prev;

	Eigen::MatrixXd RHS =
		Q -
		params_penalty.alpha * Phi_q_.transpose() *
			(dotPhi_q_ * arm_->dotq_ +
			 2 * params_penalty.xi * params_penalty.w * arm_->dotPhi_ +
			 params_penalty.w * params_penalty.w * arm_->Phi_) -
		Phi_q_.transpose() * Lambda_;

	ddot_q = A_lu_.solve(RHS);

	Lambda_ += params_penalty.alpha * arm_->Phi_;

	//	cout << "lamba: " << Lambda_.transpose() << endl;

	timelog.leave("solver_ddotq.solve");

	timelog.leave("solver_ddotq");
}

/** Integrators will call this after each time step */
void CDynamicSimulator_ALi3_Dense::post_iteration(double t) {}
