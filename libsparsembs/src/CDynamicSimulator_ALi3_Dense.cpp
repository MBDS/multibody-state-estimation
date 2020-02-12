#include <sparsembs/CAssembledRigidModel.h>
#include <sparsembs/dynamic-simulators.h>

using namespace sparsembs;
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

	m_arm->buildMassMatrix_dense(m_M);
	m_M_ldlt.compute(m_M);

	const size_t nConstraints = m_arm->m_Phi.size();
	m_Lambda.setZero(nConstraints);

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

	const size_t nDepCoords = m_arm->m_q.size();

	timelog.enter("internal_integrate");

	Eigen::VectorXd Q(nDepCoords);

	const double dt2 = dt * dt;

	const Eigen::VectorXd qp_g = -(2. / dt * m_arm->m_q + m_arm->m_dotq);
	const Eigen::VectorXd qpp_g =
		-(4. / dt2 * m_arm->m_q + 4. / dt * m_arm->m_dotq + m_arm->m_ddotq);

	m_arm->m_q += dt * m_arm->m_dotq + 0.5 * dt * dt * m_arm->m_ddotq;

	m_arm->m_dotq = (2. / dt) * m_arm->m_q + qp_g;
	m_arm->m_ddotq = (4. / dt2) * m_arm->m_q + qpp_g;

	double err = 1;
	int iter = 0;

	const double tol_dyn = 1e-6;
	const int iter_max = 20;

	m_arm->update_numeric_Phi_and_Jacobians();
	m_arm->getPhi_q_dense(m_Phi_q);

	while (err > tol_dyn && iter < iter_max)
	{
		iter++;

		// phi_0 = phi(q,l,x);
		// Get "Q" (may be dynamic)
		this->build_RHS(&Q[0] /* Q */, NULL /* we don't need "c" */);
		// Q = Qg +[0;0;0;0;-k_m*(q(5)-L_0)-c_m*qp(5)];

		Eigen::VectorXd RHS =
			0.25 * dt2 *
			(m_M * m_arm->m_ddotq +
			 m_Phi_q.transpose() * params_penalty.alpha * m_arm->m_Phi +
			 m_Phi_q.transpose() * m_Lambda - Q);
		//[K,C]=evalKC(k_m, c_m);

		// f_q = M + 0.5*dt*C+0.25*dt^2*(jac'*alpha*jac+K);
		m_A = m_M +
			  0.25 * dt2 * params_penalty.alpha * m_Phi_q.transpose() * m_Phi_q;
		m_A_lu.compute(m_A);

		const Eigen::VectorXd Aq = -m_A_lu.solve(RHS);

		m_arm->m_q += Aq;
		m_arm->m_dotq = (2. / dt) * m_arm->m_q + qp_g;
		m_arm->m_ddotq = (4. / dt2) * m_arm->m_q + qpp_g;

		// phi_0 = phi(q,l,x);
		m_arm->update_numeric_Phi_and_Jacobians();
		m_arm->getPhi_q_dense(m_Phi_q);

		m_Lambda += params_penalty.alpha * m_arm->m_Phi;
		err = Aq.norm();
	}

	// cout << "iter: " << iter << endl;

	// Proyecciones en velocidad y aceleración (faltan los términos dependientes
	// del timepo, porque en este problema no hay restricciones que dependan
	// explícitamente del tiempo).
	// qp_out = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qp);
	m_arm->m_dotq = m_A_lu.solve(m_M * m_arm->m_dotq);

	// phiqpqp_0 = phiqpqp(q, qp, l);
	m_arm->m_dotPhi_q.getAsDense(m_dotPhi_q);

	// qpp_out = f_q\((M + 0.5*dt*C + 0.25*dt^2*K)*qpp -
	// 0.25*dt^2*jac'*alpha*phiqpqp_0);
	m_arm->m_ddotq = m_A_lu.solve(
		m_M * m_arm->m_ddotq - 0.25 * dt2 * params_penalty.alpha *
								   m_Phi_q.transpose() * m_dotPhi_q *
								   m_arm->m_dotq);

	timelog.leave("internal_integrate");

	return true;
}

void CDynamicSimulator_ALi3_Dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	const size_t nDepCoords = m_arm->m_q.size();

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
	m_arm->update_numeric_Phi_and_Jacobians();

	m_arm->getPhi_q_dense(m_Phi_q);
	m_A = m_M + params_penalty.alpha * m_Phi_q.transpose() * m_Phi_q;

	m_A_lu.compute(m_A);

	// Build the RHS vector:
	// RHS = Q(q,dq) - alpha * Phi_q^t* [ \dot{Phi}_q * \dot{q} + 2 * xi * omega
	// * \dot{Phi} + omega^2 * Phi ] - Phi_q^t * \lambda
	//

	m_arm->m_dotPhi_q.getAsDense(m_dotPhi_q);

	// Solve linear system:
	// -----------------------------------
	timelog.enter("solver_ddotq.solve");

	Eigen::VectorXd ddotq_next, ddotq_prev;

	Eigen::MatrixXd RHS =
		Q -
		params_penalty.alpha * m_Phi_q.transpose() *
			(m_dotPhi_q * m_arm->m_dotq +
			 2 * params_penalty.xi * params_penalty.w * m_arm->m_dotPhi +
			 params_penalty.w * params_penalty.w * m_arm->m_Phi) -
		m_Phi_q.transpose() * m_Lambda;

	ddot_q = m_A_lu.solve(RHS);

	m_Lambda += params_penalty.alpha * m_arm->m_Phi;

	//	cout << "lamba: " << m_Lambda.transpose() << endl;

	timelog.leave("solver_ddotq.solve");

	timelog.leave("solver_ddotq");
}

/** Integrators will call this after each time step */
void CDynamicSimulator_ALi3_Dense::post_iteration(double t) {}
