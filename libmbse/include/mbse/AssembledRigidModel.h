/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include "ModelDefinition.h"

namespace mbse
{
/** This struct holds the "list of instructions" for creating an actual
 * AssembledRigidModel that models a ModelDefinition */
struct TSymbolicAssembledModel
{
	const ModelDefinition& model;  //!< My "parent" model

	/** Info on Natural Coordinate DOFs in the problem (same lenth than q_) */
	std::vector<NaturalCoordinateDOF> DOFs;

	/** Additional relative coordinates */
	std::vector<RelativeDOF> rDOFs;

	TSymbolicAssembledModel(const ModelDefinition& model_) : model(model_) {}

	void clear()
	{
		DOFs.clear();
		rDOFs.clear();
	}
};

/** The read-only "topology" of an assembled multibody model: which degrees of
 * freedom exist, how they map to points, and their initial values.
 *
 * This is completely determined by a ModelDefinition and, once built, never
 * changes. Instances are therefore safe to share (as `const`) across several
 * independent AssembledRigidModel "state" objects: each of those owns its own
 * mutable q, dq, ddq, and constraint Jacobians, but all of them can point to
 * the very same AssembledRigidModelTopology without any risk of a race, e.g.
 * when different AssembledRigidModel instances built from the same topology
 * are evaluated concurrently from different threads (as GTSAM does when
 * built with TBB). \sa AssembledRigidModel::clone()
 */
class AssembledRigidModelTopology
{
   public:
	using Ptr = std::shared_ptr<const AssembledRigidModelTopology>;

	/** Builds the topology described by a symbolic assembled model. */
	static Ptr Create(const TSymbolicAssembledModel& armi);

	const ModelDefinition& mechanism_;	//!< My "parent" model

	/** Info on each Euclidean coordinate DOF in the problem
	 * Note: DOFs_.size() + rDOFs_.size() == q_.size() */
	std::vector<NaturalCoordinateDOF> DOFs_;

	/** Info on each Relative coordinate DOF in the problem
	 * Note: DOFs_.size() + rDOFs_.size() == q_.size() */
	std::vector<RelativeDOF> rDOFs_;

	/** Reverse look-up list of Euclidean points <-> DOFs in the q vector */
	std::vector<Point2ToDOF> points2DOFs_;

	/** Maps: indices in rDOFs_ ==> indices in "q_" */
	std::vector<dof_index_t> relCoordinate2Index_;

	/** Initial values for "q", of length DOFs_.size()+rDOFs_.size() */
	Eigen::VectorXd initialQ_;

	size_t nDOFs() const { return DOFs_.size() + rDOFs_.size(); }

   private:
	explicit AssembledRigidModelTopology(const ModelDefinition& model)
		: mechanism_(model)
	{
	}
};

class AssembledRigidModel
{
	friend class ModelDefinition;  // So that class can create instances of
								   // this class.

   public:
	using Ptr = std::shared_ptr<AssembledRigidModel>;

	/** Constructor: builds a fresh mutable state (q, dq, ddq, constraint
	 * Jacobians, ...) bound to the given, immutable, and potentially-shared
	 * topology.
	 */
	explicit AssembledRigidModel(AssembledRigidModelTopology::Ptr topology);

	/** Returns an independent copy of this object, sharing the same
	 * (read-only) topology but with its own private copy of the mutable state
	 * (q, dq, ddq, Q, constraint Jacobians, ...).
	 *
	 * Since evaluating a factor mutates this state as scratch space, each
	 * factor (or, more generally, each concurrent user) must own its own
	 * clone to be safely evaluated in parallel, e.g. from different threads
	 * as GTSAM does when built with TBB.
	 */
	Ptr clone() const;

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
	void copyStateFrom(const AssembledRigidModel& o);

	/** Copies the opengl object from another instance */
	void copyOpenGLRepresentationFrom(const AssembledRigidModel& o);

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
		double E_total;	 //!< Total energy (sum of all other variables)

		double E_kin;  //!< Kinetic energy
		double E_pot;  //!< Potential energy

		TEnergyValues() : E_total(0), E_kin(0), E_pot(0) {}
	};

	/** Evaluate current energy of the system. */
	void evaluateEnergy(TEnergyValues& e) const;

	struct ComputeDependentParams
	{
		ComputeDependentParams() : maxPhiNorm(1e-13), nItersMax(10) {}

		double maxPhiNorm;
		size_t nItersMax;
	};

	struct ComputeDependentResults
	{
		ComputeDependentResults() : pos_final_phi(0), ddotq(nullptr) {}

		double pos_final_phi;  //!< Output for the final Phi(q) after refining
							   //!< positions (only valid if update_q=true).
		Eigen::VectorXd* ddotq;	 //!< Output for ddot{q}, only used if !=nullptr
								 //!< AND the input ddotz!=nullptr
	};

	/** Update dependent coordinates, velocities and accelerations from current
	 * independent ones and current state */
	void computeDependentPosVelAcc(
		const std::vector<size_t>& z_indices, bool update_q, bool update_dq,
		const ComputeDependentParams& params,
		ComputeDependentResults& out_results,
		const Eigen::VectorXd* ddotz = nullptr);

	/** Retrieves the current coordinates of a point, which may include either
	 * fixed or variable components */
	void getPointCurrentCoords(
		const size_t pt_idx, mrpt::math::TPoint2D& pt) const;

	mrpt::math::TPoint2D getPointCurrentCoords(const size_t pt_idx) const
	{
		mrpt::math::TPoint2D pt;
		getPointCurrentCoords(pt_idx, pt);
		return pt;
	}

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
	void printCoordinates(std::ostream& o = std::cout) const;

	/** Print info about all constraints */
	void printConstraints(std::ostream& o = std::cout) const;

   private:
	/** The read-only, potentially-shared topology this state was built from.
	 * \sa clone() */
	AssembledRigidModelTopology::Ptr topology_;

	/** Created upon call to getAs3DRepresentation(), this holds a list of the
	 * 3D object associated to each body in the MBS, in the same order than in
	 * parent_.bodies_[] */
	mutable std::vector<mrpt::opengl::CRenderizable::Ptr> gl_objects_;

	Eigen::Vector3d gravity_;  //!< The gravity vector (default: [0 -9.81 0])

   public:
	/** A reference to the parent MBS. Used to access the data of bodies, etc.
	 * (Aliases topology_->mechanism_.)
	 */
	const ModelDefinition& mechanism_;

	/** Info on each Euclidean coordinate DOF in the problem
	 * Note: m_DOFs.size() + m_rDOFs.size() == m_q.size()
	 * (Aliases topology_->DOFs_.) */
	const std::vector<NaturalCoordinateDOF>& DOFs_;

	/** Reverse look-up list of Euclidean points <-> DOFs in the q vector
	 * (Aliases topology_->points2DOFs_.) */
	const std::vector<Point2ToDOF>& points2DOFs_;

	/** Info on each Relative coordinate DOF in the problem
	 * Note: m_DOFs.size() + m_rDOFs.size() == m_q.size()
	 * (Aliases topology_->rDOFs_.) */
	const std::vector<RelativeDOF>& rDOFs_;

	/** Maps: indices in rDOFs_ ==> indices in "q_"
	 * (Aliases topology_->relCoordinate2Index_.) */
	const std::vector<dof_index_t>& relCoordinate2Index_;

	/** The list of all constraints (of different kinds/classes).
	 * \note This list DOES include constant-distance constraints (not like
	 * in the original list in the parent ModelDefinition)
	 */
	std::vector<ConstraintBase::Ptr> constraints_;

	/** @name State vector itself
		@{ */
	Eigen::VectorXd q_;	 //!< State vector q with all the unknowns
	Eigen::VectorXd dotq_;	//!< Velocity vector \f$ \dot{q} \f$
	/** The previously computed acceleration vector \f$ \ddot{q} \f$  */
	Eigen::VectorXd ddotq_;

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
		const Body::TRenderParams& rp) const;

	/** Animates a 3D representation of the MBS, previously built in
	 * getAs3DRepresentation() Note that there's no need to pass any 3D object
	 * as argument to this method, since smart pointers are kept internally when
	 * calling getAs3DRepresentation(). \param rp May contain colors, etc. that
	 * want to be updated into the 3D objects. \sa getAs3DRepresentation
	 */
	void update3DRepresentation(const Body::TRenderParams& rp) const;

	/** Returns the current gravity aceleration vector, used for the bodies
	 * weights (default: [0 -9.81 0]) */
	void getGravityVector(double& gx, double& gy, double& gz) const;

	/** Changes the gravity aceleration vector, used for the bodies weights
	 * (default: [0 -9.81 0]) */
	void setGravityVector(const double gx, const double gy, const double gz);

	/** @name Solvers auxiliary methods
	 *  @{ */

	/** Assemble the MBS mass matrix "M" */
	Eigen::MatrixXd buildMassMatrix_dense() const;
	std::vector<Eigen::Triplet<double>> buildMassMatrix_sparse() const;

	/** Allocates and build a sparse representation of the Mass matrix "M". The
	 * user must free the object when not needed anymore. */
	cholmod_triplet* buildMassMatrix_sparse_CHOLMOD(cholmod_common& c) const;

	/** Assemble the MBS generalized forces "Q" vector */
	void builGeneralizedForces(Eigen::VectorXd& Q) const;

	void builGeneralizedForces(double* Q) const;

	/** Call all constraint objects and command them to update their
	 * corresponding parts in the sparse Jacobians */
	void update_numeric_Phi_and_Jacobians();

	/** See constrainst realize_operating_point(). */
	void realize_operating_point() const;

	/** @} */

   private:
	mrpt::opengl::CSetOfObjects::Ptr internal_render_ground_point(
		const Point2& pt, const Body::TRenderParams& rp) const;

   public:
	// Required for aligned mem allocator (only needed in classes containing
	// fixed-size Eigen matrices)
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};	// end class AssembledRigidModel

}  // namespace mbse
