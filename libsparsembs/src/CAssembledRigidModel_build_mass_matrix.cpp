#include <sparsembs/CAssembledModelRigid.h>
#include <sparsembs/sparsembs-utils.h>

using namespace sparsembs;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

/* -------------------------------------------------------------------
				  buildMassMatrix_sparse_CHOLMOD
  Allocates and build a sparse representation of the Mass matrix "M".
	The user must free the object when not needed anymore.
-------------------------------------------------------------------*/
cholmod_triplet* CAssembledRigidModel::buildMassMatrix_sparse_CHOLMOD(
	cholmod_common& c) const
{
	timelog.enter("buildMassMatrix_sparse_CHOLMOD");
	const size_t nDOFs = m_q.size();
	const size_t nConstr = m_Phi.size();
	const size_t DIM = 2;  // 2D, 3D

	ASSERT_(nDOFs > 0)

	const size_t estimated_nnz = DIM * DIM * (2 * nDOFs + 2 * 2 * nConstr);
	const int stype = 1;  // Symmetric, stored in upper triangular only.

	// Build in triplet form:
	timelog.enter("buildMassMatrix_sparse_CHOLMOD.triplet");

	cholmod_triplet* triplet_M = cholmod_allocate_triplet(
		nDOFs, nDOFs, estimated_nnz, stype, CHOLMOD_REAL, &c);
	ASSERT_(triplet_M)

	const std::vector<CBody>& parent_bodies = m_parent.getBodies();

	// For each body:
	for (size_t i = 0; i < parent_bodies.size(); i++)
	{
		const CBody& b = parent_bodies[i];

		Matrix2d M00, M11, M01;
		b.evaluateMassMatrix(M00, M11, M01);

		const bool p0_fixed = m_parent.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = m_parent.getPointInfo(b.points[1]).fixed;

		const size_t idx_x0 =
			this->m_points2DOFs[b.points[0]]
				.dof_x;  // Will be INVALID_DOF if it's a fixed point
		const size_t idx_x1 = this->m_points2DOFs[b.points[1]].dof_x;

		if (!p0_fixed)
			insert_submatrix_in_triplet(triplet_M, idx_x0, idx_x0, M00);
		if (!p1_fixed)
			insert_submatrix_in_triplet(triplet_M, idx_x1, idx_x1, M11);
		if (!p0_fixed && !p1_fixed)
			insert_submatrix_in_triplet(
				triplet_M, std::min(idx_x0, idx_x1), std::max(idx_x0, idx_x1),
				M01);

	}  // end for each body

	timelog.leave("buildMassMatrix_sparse_CHOLMOD.triplet");

	// Convert to compressed form:
	// timelog.enter("buildMassMatrix_sparse_CHOLMOD.ccs");
	// cholmod_sparse *M = cholmod_triplet_to_sparse(triplet_M, triplet_M->nnz,
	// &c); timelog.leave("buildMassMatrix_sparse_CHOLMOD.ccs"); ASSERT_(M)
	// cholmod_free_triplet(&triplet_M, &c); // Free triplet form

	timelog.leave("buildMassMatrix_sparse_CHOLMOD");

	return triplet_M;
}

/* -------------------------------------------------------------------
				  buildMassMatrix_dense
-------------------------------------------------------------------*/
void CAssembledRigidModel::buildMassMatrix_dense(Eigen::MatrixXd& M) const
{
	timelog.enter("buildMassMatrix_dense");

	const size_t nDOFs = m_q.size();
	ASSERT_(nDOFs > 0)

	// Assemble all mass matrices:
	M.setZero(nDOFs, nDOFs);

	const std::vector<CBody>& parent_bodies = m_parent.getBodies();

	// For each body:
	for (size_t i = 0; i < parent_bodies.size(); i++)
	{
		const CBody& b = parent_bodies[i];

		const Matrix2d& M00 = b.getM00();
		const Matrix2d& M11 = b.getM11();
		const Matrix2d& M01 = b.getM01();

		const bool p0_fixed = m_parent.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = m_parent.getPointInfo(b.points[1]).fixed;

		const size_t idx_x0 =
			this->m_points2DOFs[b.points[0]]
				.dof_x;  // Will be INVALID_DOF if it's a fixed point
		const size_t idx_x1 = this->m_points2DOFs[b.points[1]].dof_x;

		if (!p0_fixed) M.block<2, 2>(idx_x0, idx_x0) += M00;
		if (!p1_fixed) M.block<2, 2>(idx_x1, idx_x1) += M11;
		if (!p0_fixed && !p1_fixed)
		{
			M.block<2, 2>(std::min(idx_x0, idx_x1), std::max(idx_x0, idx_x1)) =
				M01;
			M.block<2, 2>(std::max(idx_x0, idx_x1), std::min(idx_x0, idx_x1)) =
				M01.transpose();
		}
	}  // end for each body

	timelog.leave("buildMassMatrix_dense");
}

/* -------------------------------------------------------------------
				  buildMassMatrix_sparse
-------------------------------------------------------------------*/
void CAssembledRigidModel::buildMassMatrix_sparse(
	std::vector<Eigen::Triplet<double>>& tri) const
{
	timelog.enter("buildMassMatrix_sparse");

	const size_t nDOFs = m_q.size();
	ASSERT_(nDOFs > 0)

	// Assemble all mass matrices:
	tri.clear();

	const std::vector<CBody>& parent_bodies = m_parent.getBodies();

	// For each body:
	for (size_t i = 0; i < parent_bodies.size(); i++)
	{
		const CBody& b = parent_bodies[i];

		Matrix2d M00, M11, M01;
		b.evaluateMassMatrix(M00, M11, M01);

		const bool p0_fixed = m_parent.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = m_parent.getPointInfo(b.points[1]).fixed;

		const size_t idx_x0 =
			this->m_points2DOFs[b.points[0]]
				.dof_x;  // Will be INVALID_DOF if it's a fixed point
		const size_t idx_x1 = this->m_points2DOFs[b.points[1]].dof_x;

		if (!p0_fixed) insert_submatrix(tri, idx_x0, idx_x0, M00);
		if (!p1_fixed) insert_submatrix(tri, idx_x1, idx_x1, M11);
		if (!p0_fixed && !p1_fixed)
		{
			insert_submatrix(
				tri, std::min(idx_x0, idx_x1), std::max(idx_x0, idx_x1), M01);
			insert_submatrix(
				tri, std::max(idx_x0, idx_x1), std::min(idx_x0, idx_x1),
				M01.transpose());
		}
	}  // end for each body

	timelog.leave("buildMassMatrix_sparse");
}
