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
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense LU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Lagrange_LU_dense::CDynamicSimulator_Lagrange_LU_dense(
	const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Lagrange_LU_dense::internal_prepare()
{
	timelog.enter("solver_prepare");

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	arm_->buildMassMatrix_dense(mass_);

	timelog.leave("solver_prepare");
}

void CDynamicSimulator_Lagrange_LU_dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	timelog.enter("solver_ddotq");

	// [   M    Phi_q^t  ] [ ddot_q ] = [ Q ]
	// [ Phi_q     0     ] [ lambda ]   [ c ]
	//
	// c = - \dot{Phi_t} - \dot{Phi_q} * \dot{q}
	//  normally =>  c = - \dot{Phi_q} * \dot{q}
	//
	const size_t nDOFs = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();
	const size_t nTot = nDOFs + nConstraints;

	// Build the dense augmented matrix:
	Eigen::MatrixXd A(nTot, nTot);
	A.block(0, 0, nDOFs, nDOFs) = mass_;
	A.block(0, nDOFs, nDOFs, nConstraints).setZero();
	A.block(nDOFs, 0, nConstraints, nDOFs).setZero();
	A.block(nDOFs, nDOFs, nConstraints, nConstraints).setZero();

	// Update numeric values of the constraint Jacobians:
	timelog.enter("solver_ddotq.update_jacob");
	arm_->update_numeric_Phi_and_Jacobians();

	for (size_t i = 0; i < nConstraints; i++)
	{
		// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
		const CompressedRowSparseMatrix::row_t& row_i = arm_->Phi_q_.matrix[i];
		for (const auto& kv : row_i)
		{
			const size_t col = kv.first;
			// Insert at (col,i) because it's tranposed:

			A.coeffRef(col, nDOFs + i) = kv.second;
			A.coeffRef(nDOFs + i, col) = kv.second;
		}
	}
	timelog.leave("solver_ddotq.update_jacob");

	// Build the RHS vector:
	// --------------------------
	timelog.enter("solver_ddotq.build_rhs");
	Eigen::VectorXd RHS(nTot);
	this->build_RHS(&RHS[0], &RHS[nDOFs]);
	timelog.leave("solver_ddotq.build_rhs");

	// Solve linear system (using LU dense decomposition):
	// -------------------------------------------------------------
	timelog.enter("solver_ddotq.solve");
	const Eigen::VectorXd solution = A.partialPivLu().solve(RHS);
	timelog.leave("solver_ddotq.solve");

	ddot_q = solution.head(nDOFs);
	if (lagrangre) *lagrangre = solution.tail(nConstraints);

#if 0
	// A.saveToTextFile("A.txt");
	// RHS.saveToTextFile("RHS.txt");

	cout << "q: " << arm_->q_.transpose() << endl;
	cout << "qdot: " << arm_->dotq_.transpose() << endl;
	cout << "A:\n" << A << endl;
	cout << "RHS:\n" << RHS << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	mrpt::system::pause();
#endif

	timelog.leave("solver_ddotq");
}
