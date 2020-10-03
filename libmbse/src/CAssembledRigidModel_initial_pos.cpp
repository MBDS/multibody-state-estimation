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
#include <mbse/mbse-utils.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

/** Solves the "initial position" problem: iterates refining the position until
 * the constraints are minimized */
double CAssembledRigidModel::refinePosition(
	const double maxPhiNorm, const size_t nItersMax)
{
	timelog().enter("refinePosition");

	Eigen::MatrixXd Phi_q;
	this->update_numeric_Phi_and_Jacobians();

	size_t iter = 0;
	double phi_norm = Phi_.norm();

	timelog().registerUserMeasure("refinePosition.init_phi_norm", phi_norm);

	Eigen::FullPivLU<Eigen::MatrixXd> lu_Phiq;
	bool rebuild_lu = true;

	// Non-linear Newton iterations:
	for (; iter < nItersMax && phi_norm > maxPhiNorm; iter++)
	{
		if (rebuild_lu)
		{
			this->getPhi_q_dense(Phi_q);

			lu_Phiq.compute(Phi_q);
			rebuild_lu = false;
		}

		// Solve for increment:
		const Eigen::VectorXd q_incr = lu_Phiq.solve(Phi_);
		q_ -= q_incr;

		// Re-evaluate error:
		this->update_numeric_Phi_and_Jacobians();

		const double new_phi_norm = Phi_.norm();

		// Selective re-evaluation of the Jacobian:
		if (new_phi_norm > 1e-6) rebuild_lu = true;

		phi_norm = new_phi_norm;
	}

	timelog().registerUserMeasure("refinePosition.num_iters", iter);

	timelog().leave("refinePosition");

	return phi_norm;
}

/** Solves the "finite displacement" problem: iterates refining the position
 * until the constraints are minimized, keeping q[z_indices] fixed. */
double CAssembledRigidModel::finiteDisplacement(
	const std::vector<size_t>& z_indices, const double maxPhiNorm,
	const size_t nItersMax, bool also_correct_velocities,
	std::vector<size_t>* out_idxs_d)
{
	timelog().enter("finiteDisplacement");

	Eigen::MatrixXd Phi_q;
	this->update_numeric_Phi_and_Jacobians();

	size_t iter = 0;
	double phi_norm = Phi_.norm();

	timelog().registerUserMeasure("finiteDisplacement.init_phi_norm", phi_norm);

	std::vector<bool> q_fixed;
	q_fixed.assign(q_.size(), false);
	for (size_t i = 0; i < z_indices.size(); i++) q_fixed[z_indices[i]] = true;

	const size_t nDepCoords = q_.size() - z_indices.size();
	std::vector<size_t> idxs_d;  // make a list with the rest of indices
	idxs_d.reserve(nDepCoords);

	for (int i = 0; i < q_.size(); i++)
		if (!q_fixed[i]) idxs_d.push_back(i);

	Eigen::FullPivLU<Eigen::MatrixXd> lu_Phiq;
	bool rebuild_lu = true;

	// Non-linear Newton iterations:
	for (; iter < nItersMax && phi_norm > maxPhiNorm; iter++)
	{
		if (rebuild_lu)
		{
			this->getPhi_q_dense(Phi_q);

			mbse::removeColumns(Phi_q, z_indices);

			lu_Phiq.compute(Phi_q);
			rebuild_lu = false;
		}
		// Solve for increment:
		const Eigen::VectorXd qi_incr = lu_Phiq.solve(Phi_);

		for (size_t i = 0; i < nDepCoords; i++) q_[idxs_d[i]] -= qi_incr[i];

		// Re-evaluate error:
		this->update_numeric_Phi_and_Jacobians();

		const double new_phi_norm = Phi_.norm();

		// Selective re-evaluation of the Jacobian:
		if (new_phi_norm > 1e-6) rebuild_lu = true;

		phi_norm = new_phi_norm;
	}

	timelog().registerUserMeasure("finiteDisplacement.num_iters", iter);

	timelog().leave("finiteDisplacement");

	// Correct dependent velocities
	// --------------------------------
	if (also_correct_velocities)
	{
		timelog().enter("finiteDisplacement.dotq");

		Eigen::MatrixXd Phi_q;
		this->getPhi_q_dense(Phi_q);

		// qd = Phi_d \ (-Phi_i * dot{q}_i)
		//      -------------v-------------
		//              = vector "p"
		const size_t nConstr = Phi_q.rows();
		Eigen::VectorXd p(nConstr);
		for (size_t i = 0; i < nConstr; i++)
		{
			double r = 0;
			for (size_t j = 0; j < z_indices.size(); j++)
				r -= Phi_q(i, z_indices[j]) * dotq_[z_indices[j]];

			p[i] = r;
		}

		mbse::removeColumns(Phi_q, z_indices);

		const Eigen::VectorXd dotq_d = Phi_q.lu().solve(p);

		for (size_t i = 0; i < idxs_d.size(); i++)
			dotq_[idxs_d[i]] = dotq_d[i];

		timelog().leave("finiteDisplacement.dotq");
	}

	// Return this precomputed list of dependent indices, to save time in the
	// caller function.
	if (out_idxs_d) out_idxs_d->swap(idxs_d);

	return phi_norm;
}

/** Update dependent coordinates, velocities and accelerations from current
 * independent ones and current state */
void CAssembledRigidModel::computeDependentPosVelAcc(
	const std::vector<size_t>& z_indices, bool update_q, bool update_dq,
	const TComputeDependentParams& params,
	TComputeDependentResults& out_results, const Eigen::VectorXd* ptr_ddotz)
{
	timelog().enter("computeDependentPosVelAcc");

	// Build list of coordinates indices:
	std::vector<bool> q_fixed;
	q_fixed.assign(q_.size(), false);
	for (size_t i = 0; i < z_indices.size(); i++) q_fixed[z_indices[i]] = true;

	const size_t nDepCoords = q_.size() - z_indices.size();
	std::vector<size_t> idxs_d;  // make a list with the rest of indices
	idxs_d.reserve(nDepCoords);

	for (int i = 0; i < q_.size(); i++)
		if (!q_fixed[i]) idxs_d.push_back(i);

	// ------------------------------------------
	// Update q
	// ------------------------------------------
	if (update_q)
	{
		Eigen::MatrixXd Phi_q;
		this->update_numeric_Phi_and_Jacobians();

		size_t iter = 0;
		double phi_norm = Phi_.norm();

		timelog().registerUserMeasure(
			"computeDependentPosVelAcc.init_phi_norm", phi_norm);

		Eigen::FullPivLU<Eigen::MatrixXd> lu_Phiq;
		bool rebuild_lu = true;

		for (; iter < params.nItersMax && phi_norm > params.maxPhiNorm; iter++)
		{
			if (rebuild_lu)
			{
				this->getPhi_q_dense(Phi_q);

				mbse::removeColumns(Phi_q, z_indices);

				lu_Phiq.compute(Phi_q);
				rebuild_lu = false;
			}
			// Solve for increment:
			const Eigen::VectorXd qi_incr = lu_Phiq.solve(Phi_);

			for (size_t i = 0; i < nDepCoords; i++)
				q_[idxs_d[i]] -= qi_incr[i];

			// Re-evaluate error:
			this->update_numeric_Phi_and_Jacobians();

			const double new_phi_norm = Phi_.norm();

			// Selective re-evaluation of the Jacobian:
			if (new_phi_norm > 1e-6) rebuild_lu = true;

			phi_norm = new_phi_norm;
		}

		timelog().registerUserMeasure(
			"computeDependentPosVelAcc.num_iters", iter);
	}

	// ------------------------------------------
	// Update \dot{q}
	// ------------------------------------------
	if (update_dq)
	{
		timelog().enter("computeDependentPosVelAcc.dotq");

		Eigen::MatrixXd Phi_q;
		this->getPhi_q_dense(Phi_q);

		// qd = Phi_d \ (-Phi_i * dot{q}_i)
		//      -------------v-------------
		//              = vector "p"
		const size_t nConstr = Phi_q.rows();
		Eigen::VectorXd p(nConstr);
		for (size_t i = 0; i < nConstr; i++)
		{
			double r = 0;
			for (size_t j = 0; j < z_indices.size(); j++)
				r -= Phi_q(i, z_indices[j]) * dotq_[z_indices[j]];

			p[i] = r;
		}

		mbse::removeColumns(Phi_q, z_indices);
		const Eigen::VectorXd dotq_d = Phi_q.lu().solve(p);

		for (size_t i = 0; i < idxs_d.size(); i++)
			dotq_[idxs_d[i]] = dotq_d[i];

		timelog().leave("computeDependentPosVelAcc.dotq");
	}

	// ------------------------------------------
	// Update \ddot{q}
	// ------------------------------------------
	ASSERT_(
		(ptr_ddotz && out_results.ddotq) || (!ptr_ddotz && !out_results.ddotq));
	if (ptr_ddotz)
	{
		timelog().enter("computeDependentPosVelAcc.ddotq");

		const Eigen::VectorXd& ddotz = *ptr_ddotz;
		Eigen::VectorXd& ddotq = *out_results.ddotq;

		const size_t nConstraints = this->Phi_.size();

		Eigen::MatrixXd Phiq;
		this->getPhi_q_dense(Phiq);

		// ddot{qd} = Phiq_d \ (-Phiq_i * ddot{q}_i - dot{Phi_q} * dotq)
		//                     ----------------------v-------------------
		//                                 = vector "p"
		Eigen::VectorXd p(nConstraints);
		ASSERT_EQUAL_((int)ddotz.size(), (int)z_indices.size());

		for (size_t i = 0; i < nConstraints; i++)
		{
			double r = 0;
			// Part 1: -Phiq_i * ddot{q}_i
			for (size_t j = 0; j < z_indices.size(); j++)
				r -= Phiq(i, z_indices[j]) * ddotz[j];

			// Part 2: - dot{Phi_q} * dotq)
			const CompressedRowSparseMatrix::row_t& row_i =
				dotPhi_q_.matrix[i];
			for (CompressedRowSparseMatrix::row_t::const_iterator it =
					 row_i.begin();
				 it != row_i.end(); ++it)
			{
				const size_t col = it->first;
				r -= it->second * dotq_[col];
			}

			p[i] = r;
		}

		mbse::removeColumns(Phiq, z_indices);
		const Eigen::VectorXd ddotq_d = Phiq.lu().solve(p);

		// ------------------------------------
		// Store accelerations:
		//  ddotq[ z_indices ] <- ddotz
		//  ddotq[ idxs_d ]       <- ddotq_d
		// ------------------------------------
		ddotq.resize(q_.size());
		for (size_t i = 0; i < z_indices.size(); i++)
			ddotq[z_indices[i]] = ddotz[i];

		for (size_t i = 0; i < idxs_d.size(); i++)
			ddotq[idxs_d[i]] = ddotq_d[i];

		timelog().leave("computeDependentPosVelAcc.ddotq");
	}

	timelog().leave("computeDependentPosVelAcc");
}
