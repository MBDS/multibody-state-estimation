/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/AssembledRigidModel.h>
#include <mbse/mbse-utils.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

/* -------------------------------------------------------------------
				  buildMassMatrix_sparse_CHOLMOD
  Allocates and build a sparse representation of the Mass matrix "M".
	The user must free the object when not needed anymore.
-------------------------------------------------------------------*/
cholmod_triplet* AssembledRigidModel::buildMassMatrix_sparse_CHOLMOD(
	cholmod_common& c) const
{
	const auto tle = mrpt::system::CTimeLoggerEntry(
		timelog(), "buildMassMatrix_sparse_CHOLMOD");

	const size_t nDOFs = q_.size();
	const size_t nConstr = Phi_.size();
	const size_t DIM = 2;  // 2D, 3D

	ASSERT_(nDOFs > 0);

	const size_t estimated_nnz = DIM * DIM * (2 * nDOFs + 2 * 2 * nConstr);
	const int stype = 1;  // Symmetric, stored in upper triangular only.

	// Build in triplet form:
	timelog().enter("buildMassMatrix_sparse_CHOLMOD.triplet");

	cholmod_triplet* triplet_M = cholmod_allocate_triplet(
		nDOFs, nDOFs, estimated_nnz, stype, CHOLMOD_REAL, &c);
	ASSERT_(triplet_M);

	// For each body:
	for (const Body& b : mechanism_.bodies())
	{
		Matrix2d M00, M11, M01;
		b.evaluateMassMatrix(M00, M11, M01);

		const bool p0_fixed = mechanism_.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = mechanism_.getPointInfo(b.points[1]).fixed;

		const dof_index_t idx_x0 =
			points2DOFs_[b.points[0]]
				.dof_x;	 // Will be INVALID_DOF if it's a fixed point
		const dof_index_t idx_x1 = points2DOFs_[b.points[1]].dof_x;

		if (!p0_fixed)
			insert_submatrix_in_triplet(triplet_M, idx_x0, idx_x0, M00);
		if (!p1_fixed)
			insert_submatrix_in_triplet(triplet_M, idx_x1, idx_x1, M11);
		if (!p0_fixed && !p1_fixed)
			insert_submatrix_in_triplet(
				triplet_M, std::min(idx_x0, idx_x1), std::max(idx_x0, idx_x1),
				M01);

	}  // end for each body

	timelog().leave("buildMassMatrix_sparse_CHOLMOD.triplet");

	// Convert to compressed form:
	// timelog().enter("buildMassMatrix_sparse_CHOLMOD.ccs");
	// cholmod_sparse *M = cholmod_triplet_to_sparse(triplet_M, triplet_M->nnz,
	// &c); timelog().leave("buildMassMatrix_sparse_CHOLMOD.ccs"); ASSERT_(M)
	// cholmod_free_triplet(&triplet_M, &c); // Free triplet form

	return triplet_M;
}

/* -------------------------------------------------------------------
				  buildMassMatrix_dense
-------------------------------------------------------------------*/
Eigen::MatrixXd AssembledRigidModel::buildMassMatrix_dense() const
{
	const auto tle =
		mrpt::system::CTimeLoggerEntry(timelog(), "buildMassMatrix_dense");

	Eigen::MatrixXd M;

	const size_t nDOFs = q_.size();
	ASSERT_(nDOFs > 0);

	// Assemble all mass matrices:
	M.setZero(nDOFs, nDOFs);

	const std::vector<Body>& bodies = mechanism_.bodies();

	// For each body:
	for (const Body& b : bodies)
	{
		const Matrix2d& M00 = b.getM00();
		const Matrix2d& M11 = b.getM11();
		const Matrix2d& M01 = b.getM01();

		const bool p0_fixed = mechanism_.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = mechanism_.getPointInfo(b.points[1]).fixed;

		const dof_index_t idx_x0 =
			points2DOFs_[b.points[0]]
				.dof_x;	 // Will be INVALID_DOF if it's a fixed point
		const dof_index_t idx_x1 = points2DOFs_[b.points[1]].dof_x;

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

	return M;
}

/* -------------------------------------------------------------------
				  buildMassMatrix_sparse
-------------------------------------------------------------------*/
std::vector<Eigen::Triplet<double>>
	AssembledRigidModel::buildMassMatrix_sparse() const
{
	const auto tle =
		mrpt::system::CTimeLoggerEntry(timelog(), "buildMassMatrix_sparse");

	std::vector<Eigen::Triplet<double>> tri;

	const size_t nDOFs = q_.size();
	ASSERT_(nDOFs > 0);

	// Assemble all mass matrices:
	tri.clear();

	// For each body:
	for (const Body& b : mechanism_.bodies())
	{
		Matrix2d M00, M11, M01;
		b.evaluateMassMatrix(M00, M11, M01);

		const bool p0_fixed = mechanism_.getPointInfo(b.points[0]).fixed;
		const bool p1_fixed = mechanism_.getPointInfo(b.points[1]).fixed;

		const dof_index_t idx_x0 =
			points2DOFs_[b.points[0]]
				.dof_x;	 // Will be INVALID_DOF if it's a fixed point
		const dof_index_t idx_x1 = points2DOFs_[b.points[1]].dof_x;

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

	return tri;
}
