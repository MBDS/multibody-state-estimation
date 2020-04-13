/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include "CModelDefinition.h"

namespace mbse
{
/** This struct holds the "list of instructions" for creating an actual
 * CAssembledRigidModel that models a CModelDefinition */
struct TSymbolicAssembledModel
{
	const CModelDefinition& model;  //!< My "parent" model
	std::vector<TDOF>
		DOFs;  //!< Info on each DOF in the problem [SAME SIZE THAN m_q]

	TSymbolicAssembledModel(const CModelDefinition& model_) : model(model_) {}

	void clear() { DOFs.clear(); }
};

class CAssembledRigidModel
{
	friend class CModelDefinition;  // So that class can create instances of
									// this class.

   public:
	using Ptr = std::shared_ptr<CAssembledRigidModel>;

	/** Constructor, from a symbolic assembled model.
	 * The object "armi" can be destroyed safely after this call. The parent
	 * model in armi.model cannot.
	 */
	CAssembledRigidModel(const TSymbolicAssembledModel& armi);

	/** To be called from constraints' buildSparseStructures() methods.
	 * \return the index of the newly created row in Phi and its Jacobians.
	 */
	size_t addNewRowToConstraints();

	inline const std::vector<TPoint2DOF>& getPoints2DOFs() const
	{
		return m_points2DOFs;
	}

	/** Only to be called between objects created from the same symbolic model,
	 * this method replicates the state of "o" into "this". */
	void copyStateFrom(const CAssembledRigidModel& o);

	/** Copies the opengl object from another instance */
	void copyOpenGLRepresentationFrom(const CAssembledRigidModel& o);

	/** Solves the "initial position" problem: iterates refining the position
	 * until the constraints are minimized \return The norm of the final \Phi
	 * vector after optimization
	 */
	double refinePosition(
		const double maxPhiNorm = 1e-13, const size_t nItersMax = 10);

	/** Solves the "finite displacement" problem: iterates refining the position
	 * until the constraints are minimized, keeping q[idxs_fixed] fixed. */
	double finiteDisplacement(
		const std::vector<size_t>& idxs_fixed, const double maxPhiNorm = 1e-13,
		const size_t nItersMax = 10, bool also_correct_velocities = false,
		std::vector<size_t>* out_idxs_d = NULL);

	struct TEnergyValues
	{
		double E_total;  //!< Total energy (sum of all other variables)

		double E_kin;  //!< Kinetic energy
		double E_pot;  //!< Potential energy

		TEnergyValues() : E_total(0), E_kin(0), E_pot(0) {}
	};

	/** Evaluate current energy of the system. */
	void evaluateEnergy(TEnergyValues& e) const;

	struct TComputeDependentParams
	{
		TComputeDependentParams() : maxPhiNorm(1e-13), nItersMax(10) {}

		double maxPhiNorm;
		size_t nItersMax;
	};

	struct TComputeDependentResults
	{
		TComputeDependentResults() : pos_final_phi(0), ddotq(NULL) {}

		double pos_final_phi;  //!< Output for the final Phi(q) after refining
							   //!< positions (only valid if update_q=true).
		Eigen::VectorXd* ddotq;  //!< Output for ddot{q}, only used if !=NULL
								 //!< AND the input ddotz!=NULL
	};

	/** Update dependent coordinates, velocities and accelerations from current
	 * independent ones and current state */
	void computeDependentPosVelAcc(
		const std::vector<size_t>& z_indices, bool update_q, bool update_dq,
		const TComputeDependentParams& params,
		TComputeDependentResults& out_results,
		const Eigen::VectorXd* ddotz = NULL);

	/** Form a dense matrix from the sparse Jacobian dPhi_dq \note This method
	 * does NOT call update_numeric_Phi_and_Jacobians(), do it if required
	 * beforehand */
	inline void getPhi_q_dense(Eigen::MatrixXd& Phi_q) const
	{
		m_Phi_q.asDense(Phi_q);
	}
	inline Eigen::MatrixXd getPhi_q_dense() const
	{
		Eigen::MatrixXd m;
		m_Phi_q.asDense(m);
		return m;
	}

	/// like getPhi_q_dense() for dotPhi_q
	inline void getdotPhi_q_dense(Eigen::MatrixXd& dotPhi_q) const
	{
		m_dotPhi_q.asDense(dotPhi_q);
	}
	inline Eigen::MatrixXd getdotPhi_q_dense() const
	{
		Eigen::MatrixXd m;
		m_dotPhi_q.asDense(m);
		return m;
	}

	/// like getPhi_q_dense() for dotPhi_q
	inline void getdPhiqdq_dq_dense(Eigen::MatrixXd& dPhiqdq_dq) const
	{
		m_dPhiqdq_dq.asDense(dPhiqdq_dq);
	}
	inline Eigen::MatrixXd getdPhiqdq_dq_dense() const
	{
		Eigen::MatrixXd m;
		m_dPhiqdq_dq.asDense(m);
		return m;
	}

	/** Retrieves the current coordinates of a point, which may include either
	 * fixed or variable components */
	void getPointCurrentCoords(
		const size_t pt_idx, mrpt::math::TPoint2D& pt) const;

	/** Retrieves the current velocity of a point, which may include either
	 * fixed or variable components */
	void getPointCurrentVelocity(
		const size_t pt_idx, mrpt::math::TPoint2D& vel) const;

	/** Computes the current coordinates of a point fixed to a given body, given
	 * its relative coordinates wrt to system X:pt0->pt1, Y: orthogonal */
	void getPointOnBodyCurrentCoords(
		const size_t body_index, const mrpt::math::TPoint2D& relative_pt,
		mrpt::math::TPoint2D& out_pt) const;

   private:
	/** Created upon call to getAs3DRepresentation(), this holds a list of the
	 * 3D object associated to each body in the MBS, in the same order than in
	 * m_parent.m_bodies[] */
	mutable std::vector<mrpt::opengl::CRenderizable::Ptr> m_gl_objects;

	Eigen::Vector3d m_gravity;  //!< The gravity vector (default: [0 -9.81 0])

   public:
	/** @name State vector itself
		@{ */
	VectorXd m_q;  //!< State vector q with all the unknowns
	VectorXd m_dotq;  //!< Velocity vector \dot{q} for all the unknowns
	VectorXd m_ddotq;  //!< The previously computed acceleration vector \ddot{q}
					   //!< for all the unknowns
	/**  @} */

	/** @name Other main data, and values computed as a function of the state
	   vector.
		@{ */
	const CModelDefinition& m_parent;  //!< A reference to the parent MBS. Use
									   //!< to access the data of bodies, etc.

	std::vector<TDOF>
		m_DOFs;  //!< Info on each DOF in the problem [SAME SIZE THAN m_q]
	std::vector<TPoint2DOF> m_points2DOFs;  //!< Reverse look-up list of points
											//!< -> DOFs in the q vector

	VectorXd m_Phi;  //!< The vector of numerical values of Phi, the vector of
					 //!< constraint functions Phi=0
	VectorXd m_dotPhi;  //!< The vector of numerical values of \dot{\Phi}

	/** Jacobian dPhi_dq (as a sparse matrix) */
	TCompressedRowSparseMatrix m_Phi_q;

	/** Jacobian d(dPhi_dq)_dt (as a sparse matrix) */
	TCompressedRowSparseMatrix m_dotPhi_q;

	/** Jacobian d(Phiq*dq)_dq (as a sparse matrix) */
	TCompressedRowSparseMatrix m_dPhiqdq_dq;

	/** The list of all constraints (of different kinds/classes).
	 * \note This list DOES include constant-distance constraints (not like
	 * in the original list in the parent CModelDefinition)
	 */
	std::vector<CConstraintBase::Ptr> m_constraints;

	/** @} */

	/** Returns a 3D visualization of the model, which can be later on passed to
	 * update3DRepresentation() to animate it during a simulation. \sa
	 * update3DRepresentation \param rp May contain colors, etc. that want to be
	 * updated into the 3D objects. \note Internally, this method builds a list
	 * of 3D objects in m_gl_objects
	 */
	void getAs3DRepresentation(
		mrpt::opengl::CSetOfObjects::Ptr& outObj,
		const CBody::TRenderParams& rp) const;

	/** Animates a 3D representation of the MBS, previously built in
	 * getAs3DRepresentation() Note that there's no need to pass any 3D object
	 * as argument to this method, since smart pointers are kept internally when
	 * calling getAs3DRepresentation(). \param rp May contain colors, etc. that
	 * want to be updated into the 3D objects. \sa getAs3DRepresentation
	 */
	void update3DRepresentation(const CBody::TRenderParams& rp) const;

	/** Returns the current gravity aceleration vector, used for the bodies
	 * weights (default: [0 -9.81 0]) */
	void getGravityVector(double& gx, double& gy, double& gz) const;

	/** Changes the gravity aceleration vector, used for the bodies weights
	 * (default: [0 -9.81 0]) */
	void setGravityVector(const double gx, const double gy, const double gz);

	/** @name Solvers auxiliary methods
	 *  @{ */

	/** Assemble the MBS mass matrix "M" */
	void buildMassMatrix_dense(Eigen::MatrixXd& M) const;
	void buildMassMatrix_sparse(
		std::vector<Eigen::Triplet<double>>& M_tri) const;
	/** Allocates and build a sparse representation of the Mass matrix "M". The
	 * user must free the object when not needed anymore. */
	cholmod_triplet* buildMassMatrix_sparse_CHOLMOD(cholmod_common& c) const;

	/** Assemble the MBS generalized forces "Q" vector */
	void builGeneralizedForces(Eigen::VectorXd& Q) const;

	void builGeneralizedForces(double* Q) const;

	/** Call all constraint objects and command them to update their
	 * corresponding parts in the sparse Jacobians */
	void update_numeric_Phi_and_Jacobians();

	/** @} */

   private:
	mrpt::opengl::CSetOfObjects::Ptr internal_render_ground_point(
		const TMBSPoint& pt, const CBody::TRenderParams& rp) const;

   public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Required for aligned mem allocator (only
									 // needed in classes containing fixed-size
									 // Eigen matrices)

};  // end class CAssembledRigidModel

}  // namespace mbse
