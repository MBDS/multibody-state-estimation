#include <mbse/CAssembledRigidModel.h>
#include <mbse/dynamic-simulators.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense LU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_R_matrix_dense::CDynamicSimulator_R_matrix_dense(
	const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_R_matrix_dense::internal_prepare()
{
	timelog.enter("solver_prepare");

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	m_arm->buildMassMatrix_dense(m_mass);

	timelog.leave("solver_prepare");
}

void CDynamicSimulator_R_matrix_dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	if (lagrangre != NULL)
	{
		THROW_EXCEPTION("This class can't compute lagrange multipliers");
	}

	timelog.enter("solver_ddotq");

	// [ Phi_q ] [ ddot_q ] = [   c   ]
	// [ R^t*M ] [        ]   [ R^t*Q ]
	//
	// c = - \dot{Phi_t} - \dot{Phi_q} * \dot{q}
	//  normally =>  c = - \dot{Phi_q} * \dot{q}
	//
	const size_t nDepCoords = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();

	// Update numeric values of the constraint Jacobians:
	timelog.enter("solver_ddotq.update_jacob");

	m_arm->update_numeric_Phi_and_Jacobians();

	timelog.leave("solver_ddotq.update_jacob");

	// Get Jacobian dPhi_dq
	timelog.enter("solver_ddotq.get_dense_jacob");

	Eigen::MatrixXd Phiq(nConstraints, nDepCoords);
	m_arm->getPhi_q_dense(Phiq);

	timelog.leave("solver_ddotq.get_dense_jacob");

	// Compute R: the kernel of Phi_q
	timelog.enter("solver_ddotq.Phiq_kernel");
	Eigen::FullPivLU<Eigen::MatrixXd> lu;
	lu.compute(Phiq);
	const Eigen::MatrixXd R = lu.kernel().transpose();

	const size_t nDOFs = R.rows();

	timelog.leave("solver_ddotq.Phiq_kernel");

	ASSERT_EQUAL_(nDepCoords, nConstraints + nDOFs);

	// Build the dense augmented matrix:
	Eigen::MatrixXd A(nDepCoords, nDepCoords);
	A.block(0, 0, nConstraints, nDepCoords) = Phiq;
	A.block(nConstraints, 0, nDOFs, nDepCoords) = R;

	// Build the RHS vector:
	// --------------------------
	timelog.enter("solver_ddotq.build_rhs");
	Eigen::VectorXd RHS(nDepCoords);
	Eigen::VectorXd Q(nDepCoords);

	this->build_RHS(&Q[0], &RHS[0] /* c => [0:nConstraints-1] */);
	RHS.tail(nDOFs) = R * Q;

	timelog.leave("solver_ddotq.build_rhs");

	// Solve linear system (using LU dense decomposition):
	// -------------------------------------------------------------
	timelog.enter("solver_ddotq.solve");
	ddot_q = A.partialPivLu().solve(RHS);
	timelog.leave("solver_ddotq.solve");

#if 0
	//A.saveToTextFile("A.txt");
	//RHS.saveToTextFile("RHS.txt");

	cout << "q: " << m_arm->m_q.transpose() << endl;
	cout << "qdot: " << m_arm->m_dotq.transpose() << endl;
	cout << "A:\n" << A << endl;
	cout << "RHS:\n" << RHS << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	cout << "Phiq:\n" << Phiq << endl;
	cout << "R^t:\n" << R << endl;
	mrpt::system::pause();
#endif

	timelog.leave("solver_ddotq");
}
