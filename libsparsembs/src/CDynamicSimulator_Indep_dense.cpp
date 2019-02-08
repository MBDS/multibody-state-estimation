#include <sparsembs/CAssembledModelRigid.h>
#include <sparsembs/dynamic-simulators.h>

using namespace sparsembs;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense LU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Indep_dense::CDynamicSimulator_Indep_dense(
	const CAssembledRigidModelPtr arm_ptr)
	: CDynamicSimulatorIndepBase(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Indep_dense::internal_prepare()
{
	timelog.enter("solver_prepare");

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	m_arm->buildMassMatrix_dense(m_mass);

	timelog.leave("solver_prepare");
}

/** Compute dependent velocities and positions from the independent ones */
void CDynamicSimulator_Indep_dense::correct_dependent_q_dq()
{
	m_arm->finiteDisplacement(
		m_indep_idxs, 1e-9, 20, true /* also solve dot{q} */);
}

void CDynamicSimulator_Indep_dense::dq_plus_dz(
	const Eigen::VectorXd& dq, const Eigen::VectorXd& dz,
	Eigen::VectorXd& out_dq) const
{
	// In this model, independent accelerations are some selected indices out
	// from the vector of all dependent coordinates, so just add them:
	out_dq = dq;
	for (size_t i = 0; i < m_indep_idxs.size(); i++)
		out_dq[m_indep_idxs[i]] += dz[i];
}

// method: R matrix projection (as in section 5.2.3 of "J. García De Jalon &
// Bayo" book.
void CDynamicSimulator_Indep_dense::internal_solve_ddotz(
	double t, VectorXd& ddot_z, bool can_choose_indep_coords)
{
	timelog.enter("solver_ddotz");

	const size_t nDepCoords = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();

	// 1) Assume we know the correct values of z and \dot{z},
	//    stored in q( m_indep_idxs ) and dotq( m_indep_idxs )
	//  and find out the actual number of independent coords:
	// -----------------------------------------------------------

	// Determine number of DOFs:
	timelog.enter("solver_ddotz.update_jacob");
	m_arm->update_numeric_Phi_and_Jacobians();
	timelog.leave("solver_ddotz.update_jacob");

	// Get Jacobian dPhi_dq
	Eigen::MatrixXd Phiq(nConstraints, nDepCoords);
	timelog.enter("solver_ddotz.get_dense_jacob");
	m_arm->getPhi_q_dense(Phiq);
	timelog.leave("solver_ddotz.get_dense_jacob");

	size_t nDOFs;
	if (can_choose_indep_coords)
	{
		Eigen::FullPivLU<Eigen::MatrixXd> lu_Phiq;
		lu_Phiq.compute(Phiq);

		const size_t Phiq_rank = lu_Phiq.rank();

		// Actual # of DOFs = # dependent coords - rank of Phi_q
		nDOFs = nDepCoords - Phiq_rank;

		// Automatically choose the # of independent coordinates,
		//  by sorting the pivots of the LU decomposition (available via the
		//  permutation matrices):
		m_indep_idxs.resize(nDOFs);
		for (size_t i = 0; i < nDOFs; i++)
			m_indep_idxs[i] =
				lu_Phiq.permutationQ().indices()[nDepCoords - nDOFs + i];

		// cout << "nDOFs: " << nDOFs << ": " << m_indep_idxs << endl;
	}
	else
	{
		nDOFs = m_indep_idxs.size();
	}

	// 2) Solve for the rest of q
	// 3) Solve for rest of dotq
	// -----------------------------------------------------------
	// Assume it's already done at input!!
	// This is a task of the integrator in double
	// CDynamicSimulatorIndepBase::run(), who must assure this.

	// 4) Obtain the S & R matrices:
	//
	// [ Phi_q ]          [    b    ]
	// [ ----- ] dot{q} = [ ------  ]
	// [   B   ]          [ \dot{z} ]
	//
	//  inv(Phi_q;B) => [ S | R ]
	//
	// -----------------------------------------------------------

	Eigen::MatrixXd A(nDepCoords, nDepCoords);
	Eigen::FullPivLU<Eigen::MatrixXd> lu_A;

	A.block(0, 0, nConstraints, nDepCoords) = Phiq;
	// Fill the "B" part:
	A.block(nConstraints, 0, nDOFs, nDepCoords).setZero();
	for (size_t i = 0; i < nDOFs; i++)
		A(nConstraints + i, m_indep_idxs[i]) = 1.0;

	lu_A.compute(A);
	ASSERT_EQUAL_(lu_A.rank(), A.rows())

	const Eigen::MatrixXd A_inv = lu_A.inverse();
	const Eigen::MatrixXd S = A_inv.block(0, 0, nDepCoords, nConstraints);
	const Eigen::MatrixXd R = A_inv.block(0, nConstraints, nDepCoords, nDOFs);

	// Build the RHS vector:
	//   RHS = Rt*Q - Rt*M*Sc;
	// --------------------------
	timelog.enter("solver_ddotz.build_rhs");
	Eigen::VectorXd Q(nDepCoords);
	Eigen::VectorXd c(nConstraints);

	this->build_RHS(&Q[0], &c[0]);

	const Eigen::VectorXd RHS = R.transpose() * (Q - m_mass * S * c);

	timelog.leave("solver_ddotz.build_rhs");

	timelog.enter("solver_ddotz.solve");
	const Eigen::MatrixXd RtMR = R.transpose() * m_mass * R;
	ddot_z = RtMR.llt().solve(RHS);
	timelog.enter("solver_ddotz.solve");

#if 0
	//A.saveToTextFile("A.txt");
	//RHS.saveToTextFile("RHS.txt");

	cout << "q: " << m_arm->m_q.transpose() << endl;
	cout << "qdot: " << m_arm->m_dotq.transpose() << endl;
	cout << "A:\n" << A << endl;
	cout << "S:\n" << S << endl;
	cout << "R:\n" << R << endl;
	cout << "RHS:\n" << RHS << endl;
	cout << "solved ddotz: " << ddotz.transpose() << endl;
	cout << "solved ddot_q: " << ddot_q.transpose() << endl;
	cout << "Phiq:\n" << Phiq << endl;
	mrpt::system::pause();
#endif

	timelog.leave("solver_ddotz");
}
