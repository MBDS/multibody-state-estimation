#include <sparsembs/CAssembledRigidModel.h>
#include <sparsembs/dynamic-simulators.h>

using namespace sparsembs;
using namespace Eigen;
using namespace std;

// ---------------------------------------------------------------------------------------------
//  Solver: (Sparse) KLU
// ---------------------------------------------------------------------------------------------
CDynamicSimulator_Lagrange_KLU::CDynamicSimulator_Lagrange_KLU(
	const std::shared_ptr<CAssembledRigidModel> arm_ptr)
	: CDynamicSimulatorBase(arm_ptr),
	  ordering(orderCOLAMD),  // COLAMD is more efficient than AMD
	  m_numeric(NULL),
	  m_symbolic(NULL)
{
	klu_defaults(&m_common);
}

/** Prepare the linear systems and anything else required to really call
 * solve_ddotq() */
void CDynamicSimulator_Lagrange_KLU::internal_prepare()
{
	timelog.enter("solver_prepare");

	const size_t nDOFs = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();
	const size_t nTot = nDOFs + nConstraints;

	// Build mass matrix now and don't touch it anymore, since it's constant
	// with this formulation:
	m_arm->buildMassMatrix_sparse(m_mass_tri);

	// Start augmented matrix:
	m_A_tri = m_mass_tri;

	//  Add entries in the triplet form for the sparse Phi_q Jacobian.
	// -----------------------------------------------------------
	m_A_tri.reserve(
		m_A_tri.size() +
		2 * nConstraints *
			nDOFs);  // *IMPORTANT* Reserve mem at once to avoid reallocations,
					 // since we store pointers to places...

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
			const size_t idx0 = m_A_tri.size();

			m_A_tri.push_back(
				Eigen::Triplet<double>(itCol->first, nDOFs + i, 1.0));
			m_A_tri.push_back(
				Eigen::Triplet<double>(nDOFs + i, itCol->first, 1.0));

			m_A_tri_ptrs_Phi_q.push_back(
				const_cast<double*>(&m_A_tri[idx0].value()));
			m_A_tri_ptrs_Phi_q.push_back(
				const_cast<double*>(&m_A_tri[idx0 + 1].value()));
		}
	}

	// Analyze the pattern once:
	m_A.resize(nTot, nTot);
	m_A.setFromTriplets(m_A_tri.begin(), m_A_tri.end());

	//   int btf ;               /* use BTF pre-ordering, or not */
	//   int ordering ;          /* 0: AMD, 1: COLAMD, 2: user P and Q,
	//                            * 3: user function */
	// m_common.ordering;
	// m_common.btf;

	/* Control [UMFPACK_ORDERING] and Info [UMFPACK_ORDERING_USED] are one of:
	 */
	switch (this->ordering)
	{
		case orderAMD:
			m_common.ordering = 0;
			break;
		case orderCOLAMD:
			m_common.ordering = 1;
			break;
		default:
		    THROW_EXCEPTION("Unknown or unsupported 'ordering' value.");
	};

	m_symbolic = klu_analyze(
		m_A.rows(), m_A.outerIndexPtr(), m_A.innerIndexPtr(), &m_common);
	if (!m_symbolic)
		THROW_EXCEPTION("Error: KLU couldn't factorize the augmented matrix.");

	// m_A.toDense().saveToTextFile("A.txt");

	timelog.leave("solver_prepare");
}

CDynamicSimulator_Lagrange_KLU::~CDynamicSimulator_Lagrange_KLU()
{
	if (m_symbolic) klu_free_symbolic(&m_symbolic, &m_common);

	if (m_numeric) klu_free_numeric(&m_numeric, &m_common);
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
	const size_t nDOFs = m_arm->m_q.size();
	const size_t nConstraints = m_arm->m_Phi.size();
	const size_t nTot = nDOFs + nConstraints;

	// Update numeric values of the constraint Jacobians:
	timelog.enter("solver_ddotq.update_jacob");
	m_arm->update_numeric_Phi_and_Jacobians();
	timelog.leave("solver_ddotq.update_jacob");

	timelog.enter("solver_ddotq.update_jacob_triplets");
	// Move the updated Jacobian values to their places in the triplet form:
	{
		size_t idx = 0;
		for (size_t i = 0; i < nConstraints; i++)
		{
			// Constraint "i" goes to column "nDOFs+i" in the augmented matrix:
			const TCompressedRowSparseMatrix::row_t row_i =
				m_arm->m_Phi_q.matrix[i];
			for (TCompressedRowSparseMatrix::row_t::const_iterator itCol =
					 row_i.begin();
				 itCol != row_i.end(); ++itCol)
			{
				*m_A_tri_ptrs_Phi_q[idx++] = itCol->second;
				*m_A_tri_ptrs_Phi_q[idx++] = itCol->second;
			}
		}
	}
	timelog.leave("solver_ddotq.update_jacob_triplets");

	// Solve numeric sparse LU:
	// -----------------------------------
	timelog.enter("solver_ddotq.ccs");
	m_A.setFromTriplets(m_A_tri.begin(), m_A_tri.end());
	timelog.leave("solver_ddotq.ccs");

	timelog.enter("solver_ddotq.numeric_factor");
	if (m_numeric) klu_free_numeric(&m_numeric, &m_common);

	m_numeric = klu_factor(
		m_A.outerIndexPtr(), m_A.innerIndexPtr(), m_A.valuePtr(), m_symbolic,
		&m_common);

	if (!m_numeric)
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

	klu_solve(m_symbolic, m_numeric, m_A.cols(), 1, &RHS[0], &m_common);

	if (m_common.status != KLU_OK)
		THROW_EXCEPTION("Error: KLU couldn't solve the linear system.");

	timelog.leave("solver_ddotq.solve");

	ddot_q = RHS.head(nDOFs);
	if (lagrangre) *lagrangre = RHS.tail(nConstraints);

#if 0
	cout << "q: " << m_arm->m_q.transpose() << endl;
	cout << "qdot: " << m_arm->m_dotq.transpose() << endl;
	cout << "RHS:\n" << RHS << endl;
	cout << "solved ddotq: " << ddot_q.transpose() << endl;
	mrpt::system::pause();
#endif

	timelog.leave("solver_ddotq");
}
