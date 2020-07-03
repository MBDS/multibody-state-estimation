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
//  Solver: (Sparse) KLU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Lagrange_KLU::CDynamicSimulator_Lagrange_KLU(
	const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr),
	  ordering(orderCOLAMD),  // COLAMD is more efficient than AMD
	  numeric_(NULL),
	  symbolic_(NULL)
{
	klu_defaults(&common_);
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Lagrange_KLU::internal_prepare()
{
	timelog.enter("solver_prepare");

	const size_t nDOFs = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();
	const size_t nTot = nDOFs + nConstraints;

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	arm_->buildMassMatrix_sparse(mass_tri_);

	// Start augmented matrix:
	A_tri_ = mass_tri_;

	//  Add entries in the triplet form for the sparse Phi_q Jacobian.
	// -----------------------------------------------------------
	A_tri_.reserve(
		A_tri_.size() +
		2 * nConstraints *
			nDOFs);  // *IMPORTANT* Reserve mem at once to avoid reallocations,
					 // since we store pointers to places...

	for (size_t i = 0; i < nConstraints; i++)
	{
		// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
		const CompressedRowSparseMatrix::row_t row_i = arm_->Phi_q_.matrix[i];
		for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
				 row_i.begin();
			 itCol != row_i.end(); ++itCol)
		{
			// We have precomputed the order in which we find the numeric
			// values, just insert at their correct place:
			const size_t idx0 = A_tri_.size();

			A_tri_.push_back(
				Eigen::Triplet<double>(itCol->first, nDOFs + i, 1.0));
			A_tri_.push_back(
				Eigen::Triplet<double>(nDOFs + i, itCol->first, 1.0));

			A_tri_ptrs_Phi_q_.push_back(
				const_cast<double*>(&A_tri_[idx0].value()));
			A_tri_ptrs_Phi_q_.push_back(
				const_cast<double*>(&A_tri_[idx0 + 1].value()));
		}
	}

	// Analyze the pattern once:
	A_.resize(nTot, nTot);
	A_.setFromTriplets(A_tri_.begin(), A_tri_.end());

	//   int btf ;               /* use BTF pre-ordering, or not */
	//   int ordering ;          /* 0: AMD, 1: COLAMD, 2: user P and Q,
	//                            * 3: user function */
	// common_.ordering;
	// common_.btf;

	/* Control [UMFPACK_ORDERING] and Info [UMFPACK_ORDERING_USED] are one of:
	 */
	switch (this->ordering)
	{
		case orderAMD:
			common_.ordering = 0;
			break;
		case orderCOLAMD:
			common_.ordering = 1;
			break;
		default:
			THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};

	symbolic_ = klu_analyze(
		A_.rows(), A_.outerIndexPtr(), A_.innerIndexPtr(), &common_);
	if (!symbolic_)
		THROW_EXCEPTION("Error: KLU couldn't factorize the augmented matrix.");

	// A_.toDense().saveToTextFile("A.txt");

	timelog.leave("solver_prepare");
}

CDynamicSimulator_Lagrange_KLU::~CDynamicSimulator_Lagrange_KLU()
{
	if (symbolic_) klu_free_symbolic(&symbolic_, &common_);

	if (numeric_) klu_free_numeric(&numeric_, &common_);
}

void CDynamicSimulator_Lagrange_KLU::internal_solve_ddotq(
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

	// Update numeric values of the constraint Jacobians:
	timelog.enter("solver_ddotq.update_jacob");
	arm_->update_numeric_Phi_and_Jacobians();
	timelog.leave("solver_ddotq.update_jacob");

	timelog.enter("solver_ddotq.update_jacob_triplets");
	// Move the updated Jacobian values to their places in the triplet form:
	{
		size_t idx = 0;
		for (size_t i = 0; i < nConstraints; i++)
		{
			// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
			const CompressedRowSparseMatrix::row_t row_i =
				arm_->Phi_q_.matrix[i];
			for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
					 row_i.begin();
				 itCol != row_i.end(); ++itCol)
			{
				*A_tri_ptrs_Phi_q_[idx++] = itCol->second;
				*A_tri_ptrs_Phi_q_[idx++] = itCol->second;
			}
		}
	}
	timelog.leave("solver_ddotq.update_jacob_triplets");

	// Solve numeric sparse LU:
	// -----------------------------------
	timelog.enter("solver_ddotq.ccs");
	A_.setFromTriplets(A_tri_.begin(), A_tri_.end());
	timelog.leave("solver_ddotq.ccs");

	timelog.enter("solver_ddotq.numeric_factor");
	if (numeric_) klu_free_numeric(&numeric_, &common_);

	numeric_ = klu_factor(
		A_.outerIndexPtr(), A_.innerIndexPtr(), A_.valuePtr(), symbolic_,
		&common_);

	if (!numeric_)
		THROW_EXCEPTION(
			"Error: KLU couldn't numeric-factorize the augmented matrix.");
	timelog.leave("solver_ddotq.numeric_factor");

	// Build the RHS vector:
	// --------------------------
	timelog.enter("solver_ddotq.build_rhs");
	Eigen::VectorXd RHS(nTot);
	this->build_RHS(&RHS[0], &RHS[nDOFs]);
	timelog.leave("solver_ddotq.build_rhs");

	// Solve linear system:
	// -----------------------------------
	timelog.enter("solver_ddotq.solve");

	// Eigen::VectorXd solution(nTot);
	// KLU leaves solution in the same place than the input RHS vector:

	klu_solve(symbolic_, numeric_, A_.cols(), 1, &RHS[0], &common_);

	if (common_.status != KLU_OK)
		THROW_EXCEPTION("Error: KLU couldn't solve the linear system.");

	timelog.leave("solver_ddotq.solve");

	ddot_q = RHS.head(nDOFs);
	if (lagrangre) *lagrangre = RHS.tail(nConstraints);

#if 0
	cout << "q: " << arm_->q_.transpose() << endl;
	cout << "qdot: " << arm_->dotq_.transpose() << endl;
	cout << "RHS:\n" << RHS << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	mrpt::system::pause();
#endif

	timelog.leave("solver_ddotq");
}
