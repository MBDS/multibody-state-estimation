/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/dynamics/dynamic-simulators.h>
#include <mbse/AssembledRigidModel.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: (Sparse) CHOLMOD
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Lagrange_CHOLMOD::CDynamicSimulator_Lagrange_CHOLMOD(
	const std::shared_ptr<AssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr),
	  ordering_M(orderAMD),
	  ordering_EEt(orderAMD),
	  mass_tri_(nullptr),
	  mass_(nullptr),
	  Lm_(nullptr),
	  Phi_q_t_tri_(nullptr)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Lagrange_CHOLMOD::internal_prepare()
{
	timelog().enter("solver_prepare");

	const size_t nDOFs = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();

	// Init CHOLMOD workspace:
	cholmod_start(&cholmod_common_);

	// Since we'll factorize only definite matrices, it's OK to factor as LL'
	// instead of LDL'
	cholmod_common_.final_ll = 1;

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	mass_tri_ = arm_->buildMassMatrix_sparse_CHOLMOD(cholmod_common_);
	ASSERT_(mass_tri_ != nullptr);

	// 1) M = Lm * Lm^t
	//  Analize Mass matrix and build symbolic decomposition:
	// ---------------------------------------
	cholmod_sparse* mass_ =
		cholmod_triplet_to_sparse(mass_tri_, mass_tri_->nnz, &cholmod_common_);
	ASSERT_(mass_ != nullptr);

	cholmod_common_.nmethods = 1;
	switch (this->ordering_M)
	{
		case orderNatural:
			cholmod_common_.method[0].ordering = CHOLMOD_NATURAL;
			break;
		case orderAMD:
			cholmod_common_.method[0].ordering = CHOLMOD_AMD;
			break;
		case orderMETIS:
			cholmod_common_.method[0].ordering = CHOLMOD_METIS;
			break;
		case orderNESDIS:
			cholmod_common_.method[0].ordering = CHOLMOD_NESDIS;
			break;
		case orderCOLAMD:
			cholmod_common_.method[0].ordering = CHOLMOD_COLAMD;
			break;
		default:
			THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};
	// cholmod_common_.postorder = TRUE;

	// save_matrix( mass_, "M.txt",  &cholmod_common_);

	Lm_ = cholmod_analyze(mass_, &cholmod_common_);
	ASSERT_(Lm_ != nullptr);

	if (cholmod_common_.status != CHOLMOD_OK)
		THROW_EXCEPTION("CHOLMOD couldn't symbolic factorize M");

	// Numeric factorization:
	cholmod_factorize(mass_, Lm_, &cholmod_common_);

	// 2) Lm * E^t = Phi_q^t
	//   Build sparse Phi_q^t: (m x n)^t = (n x m)
	// For: cholmod_spsolve(CHOLMOD_L /*Lx=b*/, Lm_, Phi_q_t_ )
	// -----------------------------------------------------------
	Phi_q_t_tri_ = cholmod_allocate_triplet(
		nDOFs, nConstraints, nDOFs * nConstraints, 0 /*unsymmetric*/,
		CHOLMOD_REAL, &cholmod_common_);
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
			const size_t idx = Phi_q_t_tri_->nnz;
			static_cast<int*>(Phi_q_t_tri_->i)[idx] = itCol->first;
			static_cast<int*>(Phi_q_t_tri_->j)[idx] = i;
			ptrs_Phi_q_t_tri_.push_back(
				static_cast<double*>(Phi_q_t_tri_->x) + idx);
			Phi_q_t_tri_->nnz++;
		}
	}

	// Allocate RHS vectors:
	Q_ =
		cholmod_allocate_dense(nDOFs, 1, nDOFs, CHOLMOD_REAL, &cholmod_common_);
	ASSERT_(Q_ != nullptr);
	c_ = cholmod_allocate_dense(
		nConstraints, 1, nConstraints, CHOLMOD_REAL, &cholmod_common_);
	ASSERT_(c_ != nullptr);

	// Build the structure of E once so we can build its symbolic decomposition
	// just once now:
	// -------------------------------------------------------------------------------------------
	// Compress sparse matrix Phi_q_t:
	cholmod_sparse* Phi_q_t = cholmod_triplet_to_sparse(
		Phi_q_t_tri_, Phi_q_t_tri_->nnz, &cholmod_common_);
	ASSERTDEB_(Phi_q_t != nullptr);

	// Solve:
	//   L   *   X   = B
	//   Lm  *  E^t  = Phi_q^t
	//
	cholmod_sparse* E_t =
		cholmod_spsolve(CHOLMOD_L /*Lx=b*/, Lm_, Phi_q_t, &cholmod_common_);
	ASSERTDEB_(E_t != nullptr);
	cholmod_sparse* E = cholmod_transpose(
		E_t, 2 /* A' complex conjugate transpose */, &cholmod_common_);
	ASSERTDEB_(E != nullptr);

	//  T = E * E^t
	//  T = Lt * Lt^t
	// Analyze of "E" means analyze "E*E^t"
	// ---------------------------------------------------------------
	cholmod_common_.nmethods = 1;
	switch (this->ordering_EEt)
	{
		case orderNatural:
			cholmod_common_.method[0].ordering = CHOLMOD_NATURAL;
			break;
		case orderAMD:
			cholmod_common_.method[0].ordering = CHOLMOD_AMD;
			break;
		case orderMETIS:
			cholmod_common_.method[0].ordering = CHOLMOD_METIS;
			break;
		case orderNESDIS:
			cholmod_common_.method[0].ordering = CHOLMOD_NESDIS;
			break;
		case orderCOLAMD:
			cholmod_common_.method[0].ordering = CHOLMOD_COLAMD;
			break;
		default:
			THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};

	Lt_ = cholmod_analyze(E, &cholmod_common_);
	ASSERT_(Lt_ != nullptr);

	if (cholmod_common_.status != CHOLMOD_OK)
		THROW_EXCEPTION("CHOLMOD couldn't symbolic factorize EE'");

	cholmod_free_sparse(&E_t, &cholmod_common_);
	cholmod_free_sparse(&E, &cholmod_common_);

	timelog().leave("solver_prepare");
}

CDynamicSimulator_Lagrange_CHOLMOD::~CDynamicSimulator_Lagrange_CHOLMOD()
{
	cholmod_free_triplet(&mass_tri_, &cholmod_common_);
	cholmod_free_sparse(&mass_, &cholmod_common_);
	cholmod_free_factor(&Lm_, &cholmod_common_);
	cholmod_free_factor(&Lt_, &cholmod_common_);

	cholmod_free_triplet(&Phi_q_t_tri_, &cholmod_common_);

	cholmod_free_dense(&Q_, &cholmod_common_);
	cholmod_free_dense(&c_, &cholmod_common_);

	// Free CHOLMOD workspace:
	cholmod_finish(&cholmod_common_);
}

void CDynamicSimulator_Lagrange_CHOLMOD::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	timelog().enter("solver_ddotq");

	// [   M    Phi_q^t  ] [ ddot_q ] = [ Q ]
	// [ Phi_q     0     ] [ lambda ]   [ c ]
	//
	// c = - \dot{Phi_t} - \dot{Phi_q} * \dot{q}
	//  normally =>  c = - \dot{Phi_q} * \dot{q}
	//
	const size_t nDOFs = arm_->q_.size();
	const size_t nConstraints = arm_->Phi_.size();

	// Update numeric values of the constraint Jacobians:
	timelog().enter("solver_ddotq.update_jacob");
	arm_->update_numeric_Phi_and_Jacobians();

	// Insert Phi_q^t Jacobian in right-top block of augmented matrix:
	{
		size_t cnt = 0;
		for (size_t i = 0; i < nConstraints; i++)
		{
			// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
			const CompressedRowSparseMatrix::row_t row_i =
				arm_->Phi_q_.matrix[i];
			for (CompressedRowSparseMatrix::row_t::const_iterator itCol =
					 row_i.begin();
				 itCol != row_i.end(); ++itCol)
			{
				// We have precomputed the order in which we find the numeric
				// values, just insert at their correct place:
				*ptrs_Phi_q_t_tri_[cnt++] = itCol->second;
			}
		}
	}
	timelog().leave("solver_ddotq.update_jacob");

	// Compress sparse matrix Phi_q_t:
	timelog().enter("solver_ddotq.ccs");
	cholmod_sparse* Phi_q_t = cholmod_triplet_to_sparse(
		Phi_q_t_tri_, Phi_q_t_tri_->nnz, &cholmod_common_);
	ASSERTDEB_(Phi_q_t != nullptr);
	timelog().leave("solver_ddotq.ccs");

	// Solve:
	//   L   *   X   = B
	//   Lm  *  E^t  = Phi_q^t
	//
	timelog().enter("solver_ddotq.solve_E");
	cholmod_sparse* E_t =
		cholmod_spsolve(CHOLMOD_L /*Lx=b*/, Lm_, Phi_q_t, &cholmod_common_);
	ASSERTDEB_(E_t != nullptr);

	cholmod_sparse* E = cholmod_transpose(
		E_t, 2 /* A' complex conjugate transpose */, &cholmod_common_);
	ASSERTDEB_(E != nullptr);
	timelog().leave("solver_ddotq.solve_E");

	//  T = E * E^t
	//  T = Lt * Lt^t
	// Numeric factorization: E*E' = Lt*Lt' --> Lt=chol(E*E')
	// ---------------------------------------------------------------
	timelog().enter("solver_ddotq.numeric_factor");
	cholmod_factorize(E, Lt_, &cholmod_common_);
	timelog().leave("solver_ddotq.numeric_factor");

	// static int k=0;
	// if (!k++) mbse::save_matrix(E,"E.txt",&cholmod_common_);

	// Update the RHS vectors:
	// --------------------------
	timelog().enter("solver_ddotq.build_rhs");
	this->build_RHS(static_cast<double*>(Q_->x), static_cast<double*>(c_->x));
	timelog().leave("solver_ddotq.build_rhs");

	timelog().enter("solver_ddotq.solve");
	// Solve: Lm x2 = Q
	cholmod_dense* x2 =
		cholmod_solve(CHOLMOD_L /*Lx=b*/, Lm_, Q_, &cholmod_common_);
	ASSERTDEB_(x2 != nullptr);

	// Solve: l2 = Lt \ (E*x2-c)
	double one[2] = {1, 0}, m1[2] = {-1, 0};  // Scalars: 1 and -1
	cholmod_sdmult(
		E_t, 1 /*transpose of Et*/, one, m1, x2, c_,
		&cholmod_common_); /* c = E*x2 - c */

	cholmod_dense* l2 =
		cholmod_solve(CHOLMOD_L /*Lx=b*/, Lt_, c_, &cholmod_common_);

	// Solve: Lt^t * l = l2
	cholmod_dense* l =
		cholmod_solve(CHOLMOD_Lt /*Ltx=b*/, Lt_, l2, &cholmod_common_);

	// Solve: x = Lm^t \ (x2-E_t*l)
	cholmod_sdmult(
		E_t, 0 /*don't transpose*/, m1, one, l, x2,
		&cholmod_common_); /* x2 = x2 + E_t*l */

	cholmod_dense* x =
		cholmod_solve(CHOLMOD_Lt /*Ltx=b*/, Lm_, x2, &cholmod_common_);
	timelog().leave("solver_ddotq.solve");

	// Copy result:
	ddot_q.resize(nDOFs);
	memcpy(&ddot_q[0], x->x, sizeof(double) * nDOFs);
	if (lagrangre)
	{
		lagrangre->resize(nConstraints);
		memcpy(
			&(*lagrangre)[0], static_cast<double*>(x->x) + nDOFs,
			sizeof(double) * nConstraints);
	}

#if 0
	save_matrix_dense(Phi_q_t,"Phi_q_t.txt",&cholmod_common_);
	save_matrix_dense(cholmod_factor_to_sparse(Lm_, &cholmod_common_),"Lm.txt",&cholmod_common_);
	save_matrix_dense(E,"E.txt",&cholmod_common_);
	save_matrix_dense(cholmod_factor_to_sparse(Lt_, &cholmod_common_),"Lt.txt",&cholmod_common_);

	cout << "q: " << arm_->q_.transpose() << endl;
	cout << "qdot: " << arm_->dotq_.transpose() << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	mrpt::system::pause();
#endif

	cholmod_free_sparse(&E_t, &cholmod_common_);
	cholmod_free_sparse(&E, &cholmod_common_);

	cholmod_free_dense(&x2, &cholmod_common_);
	cholmod_free_dense(&x, &cholmod_common_);
	cholmod_free_dense(&l, &cholmod_common_);

	timelog().leave("solver_ddotq");
}
