/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/AssembledRigidModel.h>
#include <mbse/dynamics/dynamic-simulators.h>

using namespace mbse;
using namespace Eigen;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: (Sparse) KLU with the penalty formulation
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_AugmentedLagrangian_KLU::
	CDynamicSimulator_AugmentedLagrangian_KLU(
		const std::shared_ptr<AssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBasePenalty(arm_ptr), ordering(orderAMD)
{
	klu_defaults(&common_);
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_AugmentedLagrangian_KLU::internal_prepare()
{
	timelog().enter("solver_prepare");

	const size_t nDepCoords = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();

	//
	// [ M + alpha * Phi_q^t * Phi_q ] \ddot{q} = RHS
	//
	// \--------------v-------------/
	//                =A
	//
	// RHS = Q - alpha * Phi_q^t* [ \dot{Phi}_q * \dot{q} + 2 * xi * omega *
	// \dot{q} + omega^2 * Phi  ]
	//
	// Build mass matrix (constant), and set its triplet form as the beginning
	// of the total "A" matrix:
	M_tri_ = arm_->buildMassMatrix_sparse();
	A_tri_ = M_tri_;

	//  Add entries in the triplet form for the sparse Phi_q Jacobian.
	// -----------------------------------------------------------
	A_tri_.reserve(
		A_tri_.size() + nDepCoords +
		nDepCoords *
			nDepCoords);  // *IMPORTANT* Reserve mem at once to avoid
						  // reallocations, since we store pointers to places...

	PhiqtPhi_.clear();
	PhiqtPhi_.reserve(nDepCoords * nDepCoords);

	// Note: All this could be done much more efficiently if Phi_q was stored
	// in compressed column form. But since this is only computed ONCE per
	// simulation it's probably worth leave it stay...
	for (size_t i = 0; i < nDepCoords; i++)
	{
		for (size_t j = i; j < nDepCoords; j++)
		{
			// We have to evaluate the "dot product" of the columns i and j of
			// Phi_q:
			TSparseDotProduct sdp;

			for (size_t row = 0; row < nConstraints; row++)
			{
				const CompressedRowSparseMatrix::row_t& row_r =
					arm_->Phi_q_.matrix[row];

				const double *Phi_r_i = nullptr, *Phi_r_j = nullptr;

				for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
						 row_r.begin();
					 itCol != row_r.end(); ++itCol)
				{
					const size_t col = itCol->first;
					if (col > j) break;	 // We're done in this row.
					if (col != i && col != j) continue;

					if (col == i) Phi_r_i = &(itCol->second);
					if (col == j) Phi_r_j = &(itCol->second);
				}

				// Were both Phi_q[r][i] and Phi_q[r][j] != 0??
				if (Phi_r_i && Phi_r_j)
					sdp.lst_terms.push_back(
						std::pair<const double*, const double*>(
							Phi_r_i, Phi_r_j));
			}  // end for each "row"

			// Is the product != 0?
			if (!sdp.lst_terms.empty())
			{
				// Append a new triplet entry (i,j)
				A_tri_.push_back( Eigen::Triplet<double>(i,j, 1.0 /* a dummy value, it'll be updated later on by reference */ ) );
				sdp.out_ptr1 = const_cast<double*>(&(A_tri_.back().value()));

				// And also (j,i) if i!=j:
				sdp.out_ptr2 = nullptr;
				if (i != j)
				{
					A_tri_.push_back( Eigen::Triplet<double>(j,i, 1.0 /* a dummy value, it'll be updated later on by reference */ ) );
					sdp.out_ptr2 =
						const_cast<double*>(&(A_tri_.back().value()));
				}

				PhiqtPhi_.push_back(sdp);
			}

		}  // end for "j"
	}  // end for "i"

	// Analyze the pattern once:
	A_.resize(nDepCoords, nDepCoords);
	A_.setFromTriplets(A_tri_.begin(), A_tri_.end());

	M_.resize(nDepCoords, nDepCoords);
	M_.setFromTriplets(M_tri_.begin(), M_tri_.end());

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

	// Mass matrix: factorize numerically since it's constant:
	symbolic_M_ = klu_analyze(
		M_.rows(), M_.outerIndexPtr(), M_.innerIndexPtr(), &common_);
	if (!symbolic_M_)
		THROW_EXCEPTION("Error: KLU couldn't factorize the mass matrix.");
	numeric_M_ = klu_factor(
		M_.outerIndexPtr(), M_.innerIndexPtr(), M_.valuePtr(), symbolic_M_,
		&common_);

	timelog().leave("solver_prepare");
}

CDynamicSimulator_AugmentedLagrangian_KLU::
	~CDynamicSimulator_AugmentedLagrangian_KLU()
{
	if (symbolic_) klu_free_symbolic(&symbolic_, &common_);
	if (symbolic_M_) klu_free_symbolic(&symbolic_M_, &common_);

	if (numeric_) klu_free_numeric(&numeric_, &common_);
	if (numeric_M_) klu_free_numeric(&numeric_M_, &common_);
}

void CDynamicSimulator_AugmentedLagrangian_KLU::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	const size_t nDepCoords = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();

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
	this->build_RHS(&ddotq_prev[0] /* Q */, nullptr /* we don't need "c" */);

	klu_solve(symbolic_M_, numeric_M_, M_.cols(), 1, &ddotq_prev[0], &common_);

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
	timelog().enter("solver_ddotq.update_PhiqtPhiq");
	arm_->update_numeric_Phi_and_Jacobians();

	// Move the updated Jacobian values to their places in the triplet form:
	for (size_t k = 0; k < PhiqtPhi_.size(); k++)
	{
		TSparseDotProduct& sdp = PhiqtPhi_[k];

		double res = 0;
		for (size_t i = 0; i < sdp.lst_terms.size(); i++)
			res += (*sdp.lst_terms[i].first) * (*sdp.lst_terms[i].second);

		res *= params_penalty.alpha;

		if (sdp.out_ptr1) *sdp.out_ptr1 = res;
		if (sdp.out_ptr2) *sdp.out_ptr2 = res;
	}
	timelog().leave("solver_ddotq.update_PhiqtPhiq");

	// Solve numeric sparse LU:
	// -----------------------------------
	timelog().enter("solver_ddotq.ccs");
	A_.setFromTriplets(A_tri_.begin(), A_tri_.end());
	timelog().leave("solver_ddotq.ccs");

	timelog().enter("solver_ddotq.numeric_factor");
	if (numeric_) klu_free_numeric(&numeric_, &common_);

	numeric_ = klu_factor(
		A_.outerIndexPtr(), A_.innerIndexPtr(), A_.valuePtr(), symbolic_,
		&common_);

	if (!numeric_)
		THROW_EXCEPTION(
			"Error: KLU couldn't numeric-factorize the augmented matrix.");
	timelog().leave("solver_ddotq.numeric_factor");

	// Build the RHS vector:
	// RHS = M*\ddot{q}_i -  Phi_q^t* alpha * [ \dot{Phi}_q * \dot{q} + 2 * xi *
	// omega * \dot{q} + omega^2 * Phi  ]
	//                               \ ------------------------------------v
	//                               --------------------------------------/
	//                                                                    = b
	timelog().enter("solver_ddotq.build_rhs");

	// Evaluate "b":
	// b = alpha * [ \dot{Phi}_q * \dot{q} + 2 * xi * omega * \dot{Phi} +
	// omega^2 * Phi  ]
	Eigen::VectorXd b(nConstraints);

	// \dot{Phi}_q * \dot{q}
	for (size_t r = 0; r < nConstraints; r++)
	{
		const CompressedRowSparseMatrix::row_t& row_r =
			arm_->dotPhi_q_.matrix[r];
		double res = 0;
		for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
				 row_r.begin();
			 itCol != row_r.end(); ++itCol)
			res += itCol->second * arm_->dotq_[itCol->first];
		b[r] = res;
	}

	// const Eigen::VectorXd dPhiq_dq = b;

	// 2 * xi * omega * \dot{Phi}
	const double xiw2 = 2 * params_penalty.xi * params_penalty.w;
	for (size_t r = 0; r < nConstraints; r++) b[r] += xiw2 * arm_->dotPhi_[r];

	// omega^2 * Phi
	const double w2 = params_penalty.w * params_penalty.w;
	for (size_t r = 0; r < nConstraints; r++) b[r] += w2 * arm_->Phi_[r];

	// RHS2 =  alpha * Phi_q^t * b
	Eigen::VectorXd RHS2(nDepCoords);
	RHS2.setZero();
	b *= params_penalty.alpha;
	for (size_t r = 0; r < nConstraints; r++)
	{
		const CompressedRowSparseMatrix::row_t& row_r = arm_->Phi_q_.matrix[r];
		for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
				 row_r.begin();
			 itCol != row_r.end(); ++itCol)
		{
			const size_t col = itCol->first;
			RHS2[col] += itCol->second * b[r];
		}
	}

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
		// (Directly store the RHS in the in/out vector of KLU)
		ddotq_next = (M_ * ddotq_prev) - RHS2;
		klu_solve(symbolic_, numeric_, A_.cols(), 1, &ddotq_next[0], &common_);

		if (common_.status != KLU_OK)
			THROW_EXCEPTION("Error: KLU couldn't solve the linear system.");

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
