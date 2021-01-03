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

	/** Info on Natural Coordinate DOFs in the problem (same lenth than q_) */
	std::vector<NaturalCoordinateDOF> DOFs;

	/** Additional relative coordinates */
	std::vector<RelativeDOF> rDOFs;

	TSymbolicAssembledModel(const CModelDefinition& model_) : model(model_) {}

	void clear()
	{
		DOFs.clear();
		rDOFs.clear();
	}
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

	inline const std::vector<Point2ToDOF>& getPoints2DOFs() const
	{
		return points2DOFs_;
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
		std::vector<size_t>* out_idxs_d = nullptr);

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
		TComputeDependentResults() : pos_final_phi(0), ddotq(nullptr) {}

		double pos_final_phi;  //!< Output for the final Phi(q) after refining
							   //!< positions (only valid if update_q=true).
		Eigen::VectorXd* ddotq;  //!< Output for ddot{q}, only used if !=nullptr
								 //!< AND the input ddotz!=nullptr
	};

	/** Update dependent coordinates, velocities and accelerations from current
	 * independent ones and current state */
	void computeDependentPosVelAcc(
		const std::vector<size_t>& z_indices, bool update_q, bool update_dq,
		const TComputeDependentParams& params,
		TComputeDependentResults& out_results,
		const Eigen::VectorXd* ddotz = nullptr);

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

	/** Print info for each coordinate in "q" */
	void printCoordinates(std::ostream& o = std::cout);

   private:
	/** Created upon call to getAs3DRepresentation(), this holds a list of the
	 * 3D object associated to each body in the MBS, in the same order than in
	 * parent_.bodies_[] */
	mutable std::vector<mrpt::opengl::CRenderizable::Ptr> gl_objects_;

	Eigen::Vector3d gravity_;  //!< The gravity vector (default: [0 -9.81 0])

   public:
	const CModelDefinition& parent_;  //!< A reference to the parent MBS. Use
									  //!< to access the data of bodies, etc.

	/** Info on each Euclidean coordinate DOF in the problem
	 * Note: m_DOFs.size() + m_rDOFs.size() == m_q.size() */
	std::vector<NaturalCoordinateDOF> DOFs_;

	/** Reverse look-up list of Euclidean points <-> DOFs in the q vector */
	std::vector<Point2ToDOF> points2DOFs_;

	/** Info on each Relative coordinate DOF in the problem
	 * Note: m_DOFs.size() + m_rDOFs.size() == m_q.size() */
	std::vector<RelativeDOF> rDOFs_;

	/** Maps: indices in rDOFs_ ==> indices in "q_" */
	std::vector<dof_index_t> relCoordinate2Index_;

	/** The list of all constraints (of different kinds/classes).
	 * \note This list DOES include constant-distance constraints (not like
	 * in the original list in the parent CModelDefinition)
	 */
	std::vector<CConstraintBase::Ptr> constraints_;

	/** @name State vector itself
		@{ */
	Eigen::VectorXd q_;  //!< State vector q with all the unknowns
	Eigen::VectorXd dotq_;  //!< Velocity vector \dot{q} for all the unknowns
	Eigen::VectorXd ddotq_;  //!< The previously computed acceleration vector
							 //!< \ddot{q} for all the unknowns

	/** External generalized forces (gravity NOT to be included) */
	Eigen::VectorXd Q_;
	/**  @} */

	/** @name Other vectors and matrices, computed as a function of the state
	 * vector.
	 * \note You must call update_numeric_Phi_and_Jacobians() to update
	 * all these fields after updating q,dotq, ddotq.
	 *
	 *  @{ */

	/** The vector of numerical values of Phi, the vector of constraint
	 * functions Phi=0.
	 *
	 * Dimensions: m x 1
	 */
	Eigen::VectorXd Phi_;

	/** Numerical values of \dot{\Phi}
	 *
	 * Dimensions: m x 1
	 */
	Eigen::VectorXd dotPhi_;

	void defineSparseMatricesColumnCount(const size_t nDOFs)
	{
		Phi_q_.ncols = nDOFs;
		dotPhi_q_.ncols = nDOFs;
		Phiqq_times_ddq_.ncols = nDOFs;
		dotPhiqq_times_dq_.ncols = nDOFs;
	}

	void resizeConstraintCount(const size_t m)
	{
		// Add rows:
		Phi_.resize(m);
		dotPhi_.resize(m);

		// Jacobians and related matrices:
		Phi_q_.setRowCount(m);
		dotPhi_q_.setRowCount(m);
		Phiqq_times_ddq_.setRowCount(m);
		dotPhiqq_times_dq_.setRowCount(m);
	}

	/** Jacobian dPhi_dq (as a sparse matrix)
	 *
	 * Dimensions: m x n
	 */
	CompressedRowSparseMatrix Phi_q_;

	/** Jacobian d((d Phi / dq))/dt (as a sparse matrix)
	 *
	 * Dimensions: m x n
	 */
	CompressedRowSparseMatrix dotPhi_q_;

	/** Tensor-vector product: \Phiqq \ddq
	 *
	 * Dimensions: m x n
	 */
	CompressedRowSparseMatrix Phiqq_times_ddq_;

	/** Tensor-vector product: \dotPhiqq \dq
	 *
	 * Dimensions: m x n
	 */
	CompressedRowSparseMatrix dotPhiqq_times_dq_;

	/** @} */

	/** Returns a 3D visualization of the model, which can be later on passed to
	 * update3DRepresentation() to animate it during a simulation. \sa
	 * update3DRepresentation \param rp May contain colors, etc. that want to be
	 * updated into the 3D objects. \note Internally, this method builds a list
	 * of 3D objects in gl_objects_
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
		std::vector<Eigen::Triplet<double>>& tri_) const;
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
		const Point2& pt, const CBody::TRenderParams& rp) const;

   public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // Required for aligned mem allocator (only
									 // needed in classes containing fixed-size
									 // Eigen matrices)

};  // end class CAssembledRigidModel

}  // namespace mbse
