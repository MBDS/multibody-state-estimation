/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/AssembledRigidModel.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <fstream>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: Dense LU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_R_matrix_dense::CDynamicSimulator_R_matrix_dense(
	const std::shared_ptr<AssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_R_matrix_dense::internal_prepare()
{
	auto tle = mrpt::system::CTimeLoggerEntry(timelog(), "solver_prepare");

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	mass_ = arm_->buildMassMatrix_dense();
}

void CDynamicSimulator_R_matrix_dense::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	auto tle = mrpt::system::CTimeLoggerEntry(timelog(), "solver_ddotq");

	if (lagrangre != nullptr)
	{
		THROW_EXCEPTION("This class can't compute lagrange multipliers");
	}

	// [ Phi_q ] [ ddot_q ] = [   c   ]
	// [ R^t*M ] [        ]   [ R^t*Q ]
	//
	// c = - \dot{Phi_t} - \dot{Phi_q} * \dot{q}
	//  normally =>  c = - \dot{Phi_q} * \dot{q}
	//
	const size_t nDepCoords = arm_->q_.size();
	size_t nConstraints = arm_->Phi_.size();

	// Update numeric values of the constraint Jacobians:
	timelog().enter("solver_ddotq.update_jacob");

	arm_->update_numeric_Phi_and_Jacobians();

	timelog().leave("solver_ddotq.update_jacob");

	// Get Jacobian dPhi_dq
	timelog().enter("solver_ddotq.get_dense_jacob");

	// nrows=nConstraints, ncols = nDepCoords
	Eigen::MatrixXd Phiq = arm_->Phi_q_.asDense();

	timelog().leave("solver_ddotq.get_dense_jacob");

	// Compute R: the kernel of Phi_q
	timelog().enter("solver_ddotq.Phiq_kernel");
	Eigen::FullPivLU<Eigen::MatrixXd> lu;
	lu.compute(Phiq);

	lu.setThreshold(1e-8);

	const Eigen::MatrixXd R = lu.kernel();
	const size_t nDOFs = R.cols();

	// Handle duplicated constraints:
	const size_t Phi_q_rank = lu.rank();
	if (Phi_q_rank < nConstraints)
	{
		// Remove the least "relevant" constraint, according to the LU
		// elimination:
		std::vector<size_t> rowsDontPass;

		double premultiplied_threshold =
			std::abs(lu.maxPivot()) * lu.threshold();
		for (int i = 0; i < lu.nonzeroPivots(); ++i)
		{
			bool pass =
				(std::abs(lu.matrixLU().coeff(i, i)) > premultiplied_threshold);
			const int actualIdx = lu.permutationP().indices()[i];
			if (!pass) rowsDontPass.push_back(actualIdx);
		}

		unsafeRemoveRows(Phiq, rowsDontPass);
		nConstraints = Phi_q_rank;

		ASSERT_EQUAL_(static_cast<size_t>(Phiq.rows()), Phi_q_rank);
	}

	timelog().leave("solver_ddotq.Phiq_kernel");

	ASSERT_EQUAL_(nDepCoords, nConstraints + nDOFs);

	// Build the dense augmented matrix:
	Eigen::MatrixXd A(nDepCoords, nDepCoords);
	A.block(0, 0, nConstraints, nDepCoords) = Phiq;
	A.block(nConstraints, 0, nDOFs, nDepCoords) = R.transpose() * mass_;

	// Build the RHS vector:
	// --------------------------
	timelog().enter("solver_ddotq.build_rhs");
	Eigen::VectorXd RHS(nDepCoords);
	Eigen::VectorXd Q(nDepCoords);

	this->build_RHS(&Q[0], &RHS[0] /* c => [0:nConstraints-1] */);
	RHS.tail(nDOFs) = R.transpose() * Q;

	timelog().leave("solver_ddotq.build_rhs");

	// Solve linear system (using LU dense decomposition):
	// -------------------------------------------------------------
	timelog().enter("solver_ddotq.solve");
	ddot_q = A.partialPivLu().solve(RHS);
	timelog().leave("solver_ddotq.solve");
}
