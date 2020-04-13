/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/dynamic-simulators.h>
#include <mbse/CAssembledRigidModel.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: (Sparse) CHOLMOD
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Lagrange_CHOLMOD::CDynamicSimulator_Lagrange_CHOLMOD(
	const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr),
	  ordering_M(orderAMD),
	  ordering_EEt(orderAMD),
	  m_mass_tri(NULL),
	  m_mass(NULL),
	  m_Lm(NULL),
	  m_Phi_q_t_tri(NULL)
{
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Lagrange_CHOLMOD::internal_prepare()
{
	timelog.enter("solver_prepare");

	const size_t nDOFs = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();

	// Init CHOLMOD workspace:
	cholmod_start(&m_cholmod_common);

	// Since we'll factorize only definite matrices, it's OK to factor as LL'
	// instead of LDL'
	m_cholmod_common.final_ll = 1;

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	m_mass_tri = m_arm->buildMassMatrix_sparse_CHOLMOD(m_cholmod_common);
	ASSERT_(m_mass_tri != NULL);

	// 1) M = Lm * Lm^t
	//  Analize Mass matrix and build symbolic decomposition:
	// ---------------------------------------
	cholmod_sparse* m_mass = cholmod_triplet_to_sparse(
		m_mass_tri, m_mass_tri->nnz, &m_cholmod_common);
	ASSERT_(m_mass != NULL);

	m_cholmod_common.nmethods = 1;
	switch (this->ordering_M)
	{
		case orderNatural:
			m_cholmod_common.method[0].ordering = CHOLMOD_NATURAL;
			break;
		case orderAMD:
			m_cholmod_common.method[0].ordering = CHOLMOD_AMD;
			break;
		case orderMETIS:
			m_cholmod_common.method[0].ordering = CHOLMOD_METIS;
			break;
		case orderNESDIS:
			m_cholmod_common.method[0].ordering = CHOLMOD_NESDIS;
			break;
		case orderCOLAMD:
			m_cholmod_common.method[0].ordering = CHOLMOD_COLAMD;
			break;
		default:
			THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};
	// m_cholmod_common.postorder = TRUE;

	// save_matrix( m_mass, "M.txt",  &m_cholmod_common);

	m_Lm = cholmod_analyze(m_mass, &m_cholmod_common);
	ASSERT_(m_Lm != NULL);

	if (m_cholmod_common.status != CHOLMOD_OK)
		THROW_EXCEPTION("CHOLMOD couldn't symbolic factorize M");

	// Numeric factorization:
	cholmod_factorize(m_mass, m_Lm, &m_cholmod_common);

	// 2) Lm * E^t = Phi_q^t
	//   Build sparse Phi_q^t: (m x n)^t = (n x m)
	// For: cholmod_spsolve(CHOLMOD_L /*Lx=b*/, m_Lm, m_Phi_q_t )
	// -----------------------------------------------------------
	m_Phi_q_t_tri = cholmod_allocate_triplet(
		nDOFs, nConstraints, nDOFs * nConstraints, 0 /*unsymmetric*/,
		CHOLMOD_REAL, &m_cholmod_common);
	for (size_t i = 0; i < nConstraints; i++)
	{
		// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
		const TCompressedRowSparseMatrix::row_t row_i =
			m_arm->m_Phi_q.matrix[i];
		for (TCompressedRowSparseMatrix::row_t::const_iterator itCol =
				 row_i.begin();
			 itCol != row_i.end(); ++itCol)
		{
			// We have precomputed the order in which we find the numeric
			// values, just insert at their correct place:
			const size_t idx = m_Phi_q_t_tri->nnz;
			static_cast<int*>(m_Phi_q_t_tri->i)[idx] = itCol->first;
			static_cast<int*>(m_Phi_q_t_tri->j)[idx] = i;
			m_ptrs_Phi_q_t_tri.push_back(
				static_cast<double*>(m_Phi_q_t_tri->x) + idx);
			m_Phi_q_t_tri->nnz++;
		}
	}

	// Allocate RHS vectors:
	m_Q = cholmod_allocate_dense(
		nDOFs, 1, nDOFs, CHOLMOD_REAL, &m_cholmod_common);
	ASSERT_(m_Q != NULL);
	m_c = cholmod_allocate_dense(
		nConstraints, 1, nConstraints, CHOLMOD_REAL, &m_cholmod_common);
	ASSERT_(m_c != NULL);

	// Build the structure of E once so we can build its symbolic decomposition
	// just once now:
	// -------------------------------------------------------------------------------------------
	// Compress sparse matrix Phi_q_t:
	cholmod_sparse* Phi_q_t = cholmod_triplet_to_sparse(
		m_Phi_q_t_tri, m_Phi_q_t_tri->nnz, &m_cholmod_common);
	ASSERTDEB_(Phi_q_t != NULL);

	// Solve:
	//   L   *   X   = B
	//   Lm  *  E^t  = Phi_q^t
	//
	cholmod_sparse* E_t =
		cholmod_spsolve(CHOLMOD_L /*Lx=b*/, m_Lm, Phi_q_t, &m_cholmod_common);
	ASSERTDEB_(E_t != NULL);
	cholmod_sparse* E = cholmod_transpose(
		E_t, 2 /* A' complex conjugate transpose */, &m_cholmod_common);
	ASSERTDEB_(E != NULL);

	//  T = E * E^t
	//  T = Lt * Lt^t
	// Analyze of "E" means analyze "E*E^t"
	// ---------------------------------------------------------------
	m_cholmod_common.nmethods = 1;
	switch (this->ordering_EEt)
	{
		case orderNatural:
			m_cholmod_common.method[0].ordering = CHOLMOD_NATURAL;
			break;
		case orderAMD:
			m_cholmod_common.method[0].ordering = CHOLMOD_AMD;
			break;
		case orderMETIS:
			m_cholmod_common.method[0].ordering = CHOLMOD_METIS;
			break;
		case orderNESDIS:
			m_cholmod_common.method[0].ordering = CHOLMOD_NESDIS;
			break;
		case orderCOLAMD:
			m_cholmod_common.method[0].ordering = CHOLMOD_COLAMD;
			break;
		default:
			THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};

	m_Lt = cholmod_analyze(E, &m_cholmod_common);
	ASSERT_(m_Lt != NULL);

	if (m_cholmod_common.status != CHOLMOD_OK)
		THROW_EXCEPTION("CHOLMOD couldn't symbolic factorize EE'");

	cholmod_free_sparse(&E_t, &m_cholmod_common);
	cholmod_free_sparse(&E, &m_cholmod_common);

	timelog.leave("solver_prepare");
}

CDynamicSimulator_Lagrange_CHOLMOD::~CDynamicSimulator_Lagrange_CHOLMOD()
{
	cholmod_free_triplet(&m_mass_tri, &m_cholmod_common);
	cholmod_free_sparse(&m_mass, &m_cholmod_common);
	cholmod_free_factor(&m_Lm, &m_cholmod_common);
	cholmod_free_factor(&m_Lt, &m_cholmod_common);

	cholmod_free_triplet(&m_Phi_q_t_tri, &m_cholmod_common);

	cholmod_free_dense(&m_Q, &m_cholmod_common);
	cholmod_free_dense(&m_c, &m_cholmod_common);

	// Free CHOLMOD workspace:
	cholmod_finish(&m_cholmod_common);
}

void CDynamicSimulator_Lagrange_CHOLMOD::internal_solve_ddotq(
	double t, VectorXd& ddot_q, VectorXd* lagrangre)
{
	timelog.enter("solver_ddotq");

	// [   M    Phi_q^t  ] [ ddot_q ] = [ Q ]
	// [ Phi_q     0     ] [ lambda ]   [ c ]
	//
	// c = - \dot{Phi_t} - \dot{Phi_q} * \dot{q}
	//  normally =>  c = - \dot{Phi_q} * \dot{q}
	//
	const size_t nDOFs = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();

	// Update numeric values of the constraint Jacobians:
	timelog.enter("solver_ddotq.update_jacob");
	m_arm->update_numeric_Phi_and_Jacobians();

	// Insert Phi_q^t Jacobian in right-top block of augmented matrix:
	{
		size_t cnt = 0;
		for (size_t i = 0; i < nConstraints; i++)
		{
			// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
			const TCompressedRowSparseMatrix::row_t row_i =
				m_arm->m_Phi_q.matrix[i];
			for (TCompressedRowSparseMatrix::row_t::const_iterator itCol =
					 row_i.begin();
				 itCol != row_i.end(); ++itCol)
			{
				// We have precomputed the order in which we find the numeric
				// values, just insert at their correct place:
				*m_ptrs_Phi_q_t_tri[cnt++] = itCol->second;
			}
		}
	}
	timelog.leave("solver_ddotq.update_jacob");

	// Compress sparse matrix Phi_q_t:
	timelog.enter("solver_ddotq.ccs");
	cholmod_sparse* Phi_q_t = cholmod_triplet_to_sparse(
		m_Phi_q_t_tri, m_Phi_q_t_tri->nnz, &m_cholmod_common);
	ASSERTDEB_(Phi_q_t != NULL);
	timelog.leave("solver_ddotq.ccs");

	// Solve:
	//   L   *   X   = B
	//   Lm  *  E^t  = Phi_q^t
	//
	timelog.enter("solver_ddotq.solve_E");
	cholmod_sparse* E_t =
		cholmod_spsolve(CHOLMOD_L /*Lx=b*/, m_Lm, Phi_q_t, &m_cholmod_common);
	ASSERTDEB_(E_t != NULL);

	cholmod_sparse* E = cholmod_transpose(
		E_t, 2 /* A' complex conjugate transpose */, &m_cholmod_common);
	ASSERTDEB_(E != NULL);
	timelog.leave("solver_ddotq.solve_E");

	//  T = E * E^t
	//  T = Lt * Lt^t
	// Numeric factorization: E*E' = Lt*Lt' --> Lt=chol(E*E')
	// ---------------------------------------------------------------
	timelog.enter("solver_ddotq.numeric_factor");
	cholmod_factorize(E, m_Lt, &m_cholmod_common);
	timelog.leave("solver_ddotq.numeric_factor");

	// static int k=0;
	// if (!k++) mbse::save_matrix(E,"E.txt",&m_cholmod_common);

	// Update the RHS vectors:
	// --------------------------
	timelog.enter("solver_ddotq.build_rhs");
	this->build_RHS(static_cast<double*>(m_Q->x), static_cast<double*>(m_c->x));
	timelog.leave("solver_ddotq.build_rhs");

	timelog.enter("solver_ddotq.solve");
	// Solve: Lm x2 = Q
	cholmod_dense* x2 =
		cholmod_solve(CHOLMOD_L /*Lx=b*/, m_Lm, m_Q, &m_cholmod_common);
	ASSERTDEB_(x2 != NULL);

	// Solve: l2 = Lt \ (E*x2-c)
	double one[2] = {1, 0}, m1[2] = {-1, 0};  // Scalars: 1 and -1
	cholmod_sdmult(
		E_t, 1 /*transpose of Et*/, one, m1, x2, m_c,
		&m_cholmod_common); /* c = E*x2 - c */

	cholmod_dense* l2 =
		cholmod_solve(CHOLMOD_L /*Lx=b*/, m_Lt, m_c, &m_cholmod_common);

	// Solve: Lt^t * l = l2
	cholmod_dense* l =
		cholmod_solve(CHOLMOD_Lt /*Ltx=b*/, m_Lt, l2, &m_cholmod_common);

	// Solve: x = Lm^t \ (x2-E_t*l)
	cholmod_sdmult(
		E_t, 0 /*don't transpose*/, m1, one, l, x2,
		&m_cholmod_common); /* x2 = x2 + E_t*l */

	cholmod_dense* x =
		cholmod_solve(CHOLMOD_Lt /*Ltx=b*/, m_Lm, x2, &m_cholmod_common);
	timelog.leave("solver_ddotq.solve");

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
	save_matrix_dense(Phi_q_t,"Phi_q_t.txt",&m_cholmod_common);
	save_matrix_dense(cholmod_factor_to_sparse(m_Lm, &m_cholmod_common),"Lm.txt",&m_cholmod_common);
	save_matrix_dense(E,"E.txt",&m_cholmod_common);
	save_matrix_dense(cholmod_factor_to_sparse(m_Lt, &m_cholmod_common),"Lt.txt",&m_cholmod_common);

	cout << "q: " << m_arm->m_q.transpose() << endl;
	cout << "qdot: " << m_arm->m_dotq.transpose() << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	mrpt::system::pause();
#endif

	cholmod_free_sparse(&E_t, &m_cholmod_common);
	cholmod_free_sparse(&E, &m_cholmod_common);

	cholmod_free_dense(&x2, &m_cholmod_common);
	cholmod_free_dense(&x, &m_cholmod_common);
	cholmod_free_dense(&l, &m_cholmod_common);

	timelog.leave("solver_ddotq");
}
