#include <sparsembs/CAssembledModelRigid.h>
#include <sparsembs/dynamic-simulators.h>

using namespace sparsembs;
using namespace Eigen;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense solver with the Augmented Lagrange formulation (ALF)
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_AugmentedLagrangian_Dense::
	CDynamicSimulator_AugmentedLagrangian_Dense(
		const CAssembledRigidModelPtr arm_ptr)
	: CDynamicSimulatorBasePenalty(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_AugmentedLagrangian_Dense::internal_prepare()
{
	timelog.enter("solver_prepare");

	m_arm->buildMassMatrix_dense(m_M);
	m_M_ldlt.compute(m_M);

	timelog.leave("solver_prepare");
}

CDynamicSimulator_AugmentedLagrangian_Dense::
	~CDynamicSimulator_AugmentedLagrangian_Dense()
{
}

void CDynamicSimulator_AugmentedLagrangian_Dense::internal_solve_ddotq(
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
	Eigen::VectorXd ddotq_prev(nDepCoords), ddotq_next(nDepCoords);
	this->build_RHS(&ddotq_prev[0] /* Q */, NULL /* we don't need "c" */);

	ddotq_prev = m_M_ldlt.solve(ddotq_prev);

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
	m_arm->update_numeric_Phi_and_Jacobians();

	m_arm->getPhi_q_dense(m_Phi_q);
	m_A = m_M + params_penalty.alpha * m_Phi_q.transpose() * m_Phi_q;

	m_A_lu.compute(m_A);

	// Build the RHS vector:
	// RHS = M*\ddot{q}_i -  Phi_q^t* alpha * [ \dot{Phi}_q * \dot{q} + 2 * xi *
	// omega * \dot{Phi} + omega^2 * Phi  ]
	//                               \ ------------------------------------v
	//                               --------------------------------------/
	//                                                                    = b
	timelog.enter("solver_ddotq.build_rhs");

	m_arm->m_dotPhi_q.getAsDense(m_dotPhi_q);

	const Eigen::MatrixXd RHS2 =
		params_penalty.alpha * m_Phi_q.transpose() *
		(m_dotPhi_q * m_arm->m_dotq +
		 2 * params_penalty.xi * params_penalty.w * m_arm->m_dotPhi +
		 params_penalty.w * params_penalty.w * m_arm->m_Phi);

	timelog.leave("solver_ddotq.build_rhs");

	// Solve linear system:
	// -----------------------------------
	timelog.enter("solver_ddotq.solve");

	Eigen::VectorXd RHS(nDepCoords);

	const double MAX_DDOT_INCR_NORM = 1e-4 * nDepCoords;
	const size_t MAX_ITERS = 10;

	double ddot_incr_norm;
	size_t iter = 0;
	do
	{
		// RHS = M*\ddot{q}_i - RHS2
		RHS = (m_M * ddotq_prev) - RHS2;
		ddotq_next = m_A_lu.solve(RHS);

		ddot_incr_norm = (ddotq_next - ddotq_prev).norm();
		// cout << "iter: " << iter<< endl << "prev: " << ddotq_prev.transpose()
		// << "\nnext: " << ddotq_next.transpose() << "\n  norm: " <<
		// ddot_incr_norm << endl << endl;

		ddotq_prev = ddotq_next;
	} while (ddot_incr_norm > MAX_DDOT_INCR_NORM && ++iter < MAX_ITERS);

	ddot_q.swap(ddotq_next);

	timelog.leave("solver_ddotq.solve");

	ASSERTDEBMSG_(
		((RHS.array() == RHS.array()).all()), "NaN found in result ddotq")

	timelog.leave("solver_ddotq");
}

/** Integrators will call this after each time step */
void CDynamicSimulator_AugmentedLagrangian_Dense::post_iteration(double t)
{
#if 0

	const size_t nConstraints = m_arm->m_Phi.size();

	// ---------------------------
	// Projection in q
	// ---------------------------
	const Eigen::VectorXd q0 = m_arm->m_q;
	Eigen::VectorXd Lambda(nConstraints);
	Lambda.setZero();

	Eigen::VectorXd rhs, Aq, Aqp;

	for (int i=0;i<3;i++)
	{
		// Update numeric values of the constraint Jacobians:
		m_arm->update_numeric_Phi_and_Jacobians();
		m_arm->getPhi_q_dense(m_Phi_q);

		Lambda += params_penalty.alpha * m_arm->m_Phi;

		// Solve for increments Aq:
		//
		// -> Lambda = Lambda + alpha * Phi
		// -> [M+alpha * Phi_q^t * Phi_q] Aq = -[ M (qi-q0) + Phi_q^t * Lambda ]
		rhs = m_M *( q0 - m_arm->m_q ) - m_Phi_q.transpose() * Lambda;

		Aq = m_A_lu.solve(rhs);
		m_arm->m_q += Aq;

		cout << "iter: " << i << " |Aq|=" << Aq.norm() << endl;
	}
	cout << "\n";

	// ---------------------------
	// Projection in dot{q}
	// ---------------------------
//	Eigen::FullPivLU<Eigen::MatrixXd> lu;
//	lu.compute(m_Phi_q);
//	const Eigen::MatrixXd V = lu.kernel();
//	m_arm->m_dotq = (V*V.transpose()) * m_arm->m_dotq;

	timelog.leave("solver.post_iteration");
#endif  // 0	timelog.enter("solver.post_iteration");
}
