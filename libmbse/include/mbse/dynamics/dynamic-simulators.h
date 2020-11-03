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

#include <mbse/mbse-common.h>
#include <list>

namespace mbse
{
struct TPointState
{
	TPointState(
		const mrpt::math::TPoint2D& _pos, const mrpt::math::TPoint2D _vel)
		: pos(_pos), vel(_vel)
	{
	}

	mrpt::math::TPoint2D pos, vel;
};

typedef std::pair<double, TPointState> timestamped_point_t;

enum TOrderingMethods
{
	orderNatural = 0,  //!< Leave variables in their natural order
	orderAMD,  //!< Use AMD algorithm
	orderCOLAMD,  //!< Use COLAMD algorithm (Only for KLU (??) )
	orderMETIS,  //!< Use METIS graph-based partitioning
	orderCHOLMOD,  //!< Use AMD/COLAMD then METIS
	orderNESDIS,  //!< CHOLMOD's version of nested graph disection
	orderTryKeepBest  //!< Try different methods and keep the best one
};

/** Logging structure for CDynamicSimulatorBase's "sensors" */
struct TSensorData
{
	size_t pnt_index;  //!< In the original MBS model
	const double* pos[3];  //!< Pointers to the up-to-date coordinates (X,Y,Z)
	const double* vel[3];  //!< Pointers to the up-to-date velocities (X,Y,Z)
	/** Log of sensed data: */
	std::deque<timestamped_point_t> log;
};

class CAssembledRigidModel;  //!< A MBS preprocessed and ready for
							 //!< kinematic/dynamic simulations.

enum ODE_integrator_t
{
	ODE_Euler = 0,  //!< Simple, explicit, Euler method
	ODE_Trapezoidal,  //!< Implicit 2nd order method
	ODE_RK4  //!< Explicit Runge-Kutta 4th order method
};

/** State of the simulation, passed to a user-provided function */
struct TSimulationState
{
	double t;
	const CAssembledRigidModel* const
		arm;  //!< From this object you can retrieve the current "q"
			  //!< coordinates, velocities, bodies, etc.

	TSimulationState(const CAssembledRigidModel* arm_);
};

using TSimulationStateRef = const TSimulationState&;

using simul_callback_t = std::function<void(TSimulationStateRef)>;

class CDynamicSimulatorBase;

/** The common part of all simulators */
class CDynamicSimulatorBase
{
   public:
	using Ptr = std::shared_ptr<CDynamicSimulatorBase>;

	/** A class factory, creates a dynamic simulator from a string with the
	 * class name: "CDynamicSimulator_Lagrange_LU_dense",
	 * "CDynamicSimulator_Lagrange_UMFPACK", ...
	 */
	static Ptr Create(
		const std::string& name,
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);

	CDynamicSimulatorBase(std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulatorBase();

	struct TParameters
	{
		TParameters() = default;

		/**  Method for numerical integration of ODE system */
		ODE_integrator_t ode_solver = ODE_Euler;

		/** For fixed-time integrators, the fixed time step */
		double time_step = 1e-3;

		/** Called AFTER each new simulation step */
		simul_callback_t user_callback;
	};

	TParameters params;  //!< The simulator parameters

	/** One-time preparation of the linear systems and anything else required,
	 * before starting to call solve_ddotq()
	 *  ** MUST BE CALLED BEFORE solve_ddotq() **
	 */
	void prepare();

	/** Solve for the current accelerations
	 *  You MUST call prepare() before this method.
	 */
	virtual void solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	/** Integrators will call this before solve_ddotq() once per time step */
	virtual void pre_iteration(double t) {}

	/** Integrators will call this after each time step */
	virtual void post_iteration(double t) {}

	/** Runs a dynamic simulation for a given time span
	 * \return The actual final time of the simulation, which will always be >=
	 * t_end. It may be >t_end due the fixed time step of some integrators.
	 */
	virtual double run(const double t_ini, const double t_end);

	const std::shared_ptr<CAssembledRigidModel>& get_model() const
	{
		return arm_ptr_;
	}

	CAssembledRigidModel* get_model_non_const() const { return arm_ptr_.get(); }

	/** \name Sensors
		 @{ */

	/** Add a "sensor" that grabs the position of a given point.
	 * \sa saveSensorLogsToFile
	 */
	void addPointSensor(const size_t pnt_index);

	/** Save all logged data to a text file, loadable from MATLAB with "load()".
	 * \return false on any error, true if all go ok.
	 */
	bool saveSensorLogsToFile(const std::string& filename) const;

	/** @} */

   protected:
	//!< The smart pointer. Normally use arm_ which is faster
	const std::shared_ptr<CAssembledRigidModel> arm_ptr_;
	CAssembledRigidModel* const arm_;  //!< The prepared MBS model.

	/** Build the two subvectors of the RHS of the motion equation, with the
	 * generalized forces "Q" and the Jacobian-based term "c". "Q" must have
	 * space for nDOFs doubles, and "c" must have space for "nConstraints".
	 */
	void build_RHS(double* Q, double* c);

	/** Prepare the linear systems and anything else required to really call
	 * solve_ddotq() */
	virtual void internal_prepare() = 0;

	/** Solve for the current accelerations */
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q,
		Eigen::VectorXd* lagrangre = NULL) = 0;

	/** Implement a especific combination of dynamic formulation + integrator.
	 *  \return false if it's not implemented, so it should fallback to generic
	 * integrator + internal_solve_ddotq()
	 */
	virtual bool internal_integrate(
		double t, double dt, const ODE_integrator_t integr)
	{
		return false;
	}

	// Auxiliary variables of the ODE integrators (declared here to avoid
	// reallocating mem)
	Eigen::VectorXd q0;  // Backup of state.
	Eigen::VectorXd k1, k2, k3, k4;
	Eigen::VectorXd v1, v2, v3, v4;  // \dot{q}
   private:
	Eigen::VectorXd ddotq1, ddotq2, ddotq3, ddotq4;  // \ddot{q}

   protected:
	bool init_;  //!< Used to indicate if user has called prepare()

	/** List of "sensed point_index" -> list of logged data.
	 * Updated by addPointSensor()
	 */
	std::list<TSensorData> sensors_;
};

class CDynamicSimulatorIndepBase;

/** Especialization of simulator for formulations in independent coordinates (it
 * requires different integrators) */
class CDynamicSimulatorIndepBase : public CDynamicSimulatorBase
{
   public:
	using Ptr = std::shared_ptr<CDynamicSimulatorIndepBase>;

	CDynamicSimulatorIndepBase(std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulatorIndepBase();

	/** Runs a dynamic simulation for a given time span
	 * \return The actual final time of the simulation, which will always be >=
	 * t_end. It may be >t_end due the fixed time step of some integrators.
	 */
	virtual double run(const double t_ini, const double t_end);

	/** Solve for the current independent accelerations
	 *  You MUST call prepare() before this method.
	 */
	void solve_ddotz(
		double t, Eigen::VectorXd& ddot_z, bool can_choose_indep_coords);

	/** Performs the addition of velocities: out_dq = dq +
	 * independent2dependent(dz) */
	virtual void dq_plus_dz(
		const Eigen::VectorXd& dq, const Eigen::VectorXd& dz,
		Eigen::VectorXd& out_dq) const = 0;

	/** Compute dependent velocities and positions from the independent ones */
	virtual void correct_dependent_q_dq() = 0;

	virtual const std::vector<size_t>& independent_coordinate_indices()
		const = 0;
	virtual void independent_coordinate_indices(
		const std::vector<size_t>& idxs) = 0;

   protected:
	/** Wrapper for ddotq computation, from ddotz */
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	/** Solve for the current accelerations of independent coords */
	virtual void internal_solve_ddotz(
		double t, Eigen::VectorXd& ddot_z, bool can_choose_indep_coords) = 0;

	// Auxiliary variables of the ODE integrators (declared here to avoid
	// reallocating mem)
	Eigen::VectorXd ddotz1, ddotz2, ddotz3, ddotz4;  // \ddot{z}
};

class CDynamicSimulator_Lagrange_LU_dense : public CDynamicSimulatorBase
{
   public:
	CDynamicSimulator_Lagrange_LU_dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	Eigen::MatrixXd mass_;  //!< The MBS constant mass matrix
};

class CDynamicSimulator_R_matrix_dense : public CDynamicSimulatorBase
{
   public:
	CDynamicSimulator_R_matrix_dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	Eigen::MatrixXd mass_;  //!< The MBS constant mass matrix
};

/** R matrix projection method (as in section 5.2.3 of "J. García De Jalon &
 * Bayo" book) */
class CDynamicSimulator_Indep_dense : public CDynamicSimulatorIndepBase
{
   public:
	CDynamicSimulator_Indep_dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);

	void dq_plus_dz(
		const Eigen::VectorXd& dq, const Eigen::VectorXd& dz,
		Eigen::VectorXd& out_dq) const override;
	/** Compute dependent velocities and positions from the independent ones */
	void correct_dependent_q_dq() override;

	const std::vector<size_t>& independent_coordinate_indices() const override
	{
		return indep_idxs_;
	}
	void independent_coordinate_indices(
		const std::vector<size_t>& idxs) override
	{
		indep_idxs_ = idxs;
	}

   private:
	void internal_prepare() override;
	void internal_solve_ddotz(
		double t, Eigen::VectorXd& ddot_z,
		bool can_choose_indep_coords) override;

	Eigen::MatrixXd mass_;  //!< The MBS constant mass matrix
	/** The indices in "q" of those coordinates to be used as "independent" (z)
	 */
	std::vector<size_t> indep_idxs_;
};

class CDynamicSimulator_Lagrange_CHOLMOD : public CDynamicSimulatorBase
{
   public:
	CDynamicSimulator_Lagrange_CHOLMOD(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_Lagrange_CHOLMOD();

	TOrderingMethods ordering_M;  //!< The ordering algorithm for factorizing M
	TOrderingMethods
		ordering_EEt;  //!< The ordering algorithm for factorizing E*E'

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	cholmod_common cholmod_common_;
	cholmod_triplet* mass_tri_;
	cholmod_sparse* mass_;
	cholmod_factor* Lm_;  //!< Mass = Lm * Lm'
	cholmod_factor* Lt_;  //!< E*E' = Lt*Lt'
	cholmod_triplet* Phi_q_t_tri_;
	std::vector<double*>
		ptrs_Phi_q_t_tri_;  //!< Pointers to elements in Phi_q_t_tri_ in the
							//!< same order than in the sparse Jacobian.
	cholmod_dense *Q_, *c_, *z_;  //!< RHS & auxiliary vectors
};

class CDynamicSimulator_Lagrange_UMFPACK : public CDynamicSimulatorBase
{
   public:
	CDynamicSimulator_Lagrange_UMFPACK(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_Lagrange_UMFPACK();

	TOrderingMethods ordering;

   private:
	void internal_prepare() override;
	void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q,
		Eigen::VectorXd* lagrangre = NULL) override;

	std::vector<Eigen::Triplet<double>> mass_tri_;
	std::vector<Eigen::Triplet<double>>
		A_tri_;  //!< Augmented matrix (triplet form)
	std::vector<double*> A_tri_ptrs_Phi_q_;  //!< Placeholders for the non-zero
											 //!< entries of the Jacobian.
	Eigen::SparseMatrix<double> A_;  //!< Augmented matrix (CCS)

	void* numeric_;
	void* symbolic_;

	double umf_control_[UMFPACK_CONTROL];
	double umf_info_[UMFPACK_INFO];
};

class CDynamicSimulator_Lagrange_KLU : public CDynamicSimulatorBase
{
   public:
	CDynamicSimulator_Lagrange_KLU(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_Lagrange_KLU();

	TOrderingMethods ordering;

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	std::vector<Eigen::Triplet<double>> mass_tri_;
	std::vector<Eigen::Triplet<double>>
		A_tri_;  //!< Augmented matrix (triplet form)
	std::vector<double*> A_tri_ptrs_Phi_q_;  //!< Placeholders for the non-zero
											 //!< entries of the Jacobian.
	Eigen::SparseMatrix<double> A_;  //!< Augmented matrix (CCS)

	klu_common common_;
	klu_numeric* numeric_;
	klu_symbolic* symbolic_;
};

class CDynamicSimulatorBasePenalty : public CDynamicSimulatorBase
{
   public:
	struct TPenaltyParams
	{
		double alpha;  //!< Penalty multiplier (default: 1e4)
		double w;  //!< Penalty omega (default: 10)
		double xi;  //!< Penalty xi (default: 1)

		TPenaltyParams() : alpha(1e7), w(10), xi(1) {}
	};

	TPenaltyParams
		params_penalty;  //!< Parameters of the "penalty" dynamic formulation

	CDynamicSimulatorBasePenalty(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr)
		: CDynamicSimulatorBase(arm_ptr)
	{
	}
	virtual ~CDynamicSimulatorBasePenalty() {}
};

class CDynamicSimulator_AugmentedLagrangian_KLU
	: public CDynamicSimulatorBasePenalty
{
   public:
	CDynamicSimulator_AugmentedLagrangian_KLU(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_AugmentedLagrangian_KLU();

	TOrderingMethods ordering;

	const Eigen::SparseMatrix<double>& getA() const { return A_; }

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	struct TSparseDotProduct
	{
		std::vector<std::pair<const double*, const double*>> lst_terms;
		double *out_ptr1,
			*out_ptr2;  //!< Store the result of the dot product in these
						//!< pointers, if they are not NULL.
	};

	std::vector<Eigen::Triplet<double>> A_tri_,
		M_tri_;  //!< Augmented matrix (triplet form)
	std::vector<TSparseDotProduct>
		PhiqtPhi_;  //!< Quick list of operations needed to update the product
					//!< Phi_q^t * Phi_q and store it into A_tri_.
	Eigen::SparseMatrix<double> A_, M_;  //!< Augmented matrix (CCS)

	klu_common common_;
	klu_numeric *numeric_ = nullptr, *numeric_M_ = nullptr;
	klu_symbolic *symbolic_ = nullptr, *symbolic_M_ = nullptr;
};

class CDynamicSimulator_AugmentedLagrangian_Dense
	: public CDynamicSimulatorBasePenalty
{
   public:
	CDynamicSimulator_AugmentedLagrangian_Dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_AugmentedLagrangian_Dense();

	/** Integrators will call this after each time step */
	virtual void post_iteration(double t);

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	Eigen::MatrixXd M_;  //!< The MBS constant mass matrix
	Eigen::LDLT<Eigen::MatrixXd> M_ldlt_;

	// Data updated during solve(), then reused during post_iteration():
	Eigen::MatrixXd A_, Phi_q_, dotPhi_q_;
	Eigen::FullPivLU<Eigen::MatrixXd> A_lu_;
};

class CDynamicSimulator_ALi3_Dense : public CDynamicSimulatorBasePenalty
{
   public:
	CDynamicSimulator_ALi3_Dense(
		const std::shared_ptr<CAssembledRigidModel> arm_ptr);
	virtual ~CDynamicSimulator_ALi3_Dense();

	/** Integrators will call this after each time step */
	virtual void post_iteration(double t);

   private:
	virtual void internal_prepare();
	virtual void internal_solve_ddotq(
		double t, Eigen::VectorXd& ddot_q, Eigen::VectorXd* lagrangre = NULL);

	/** Implement a especific combination of dynamic formulation + integrator.
	 *  \return false if it's not implemented, so it should fallback to generic
	 * integrator + internal_solve_ddotq()
	 */
	virtual bool internal_integrate(
		double t, double dt, const ODE_integrator_t integr);

	Eigen::MatrixXd M_;  //!< The MBS constant mass matrix
	Eigen::LDLT<Eigen::MatrixXd> M_ldlt_;

	// Data updated during solve(), then reused during post_iteration():
	Eigen::MatrixXd A_, Phi_q_, dotPhi_q_;
	Eigen::FullPivLU<Eigen::MatrixXd> A_lu_;

	Eigen::VectorXd Lambda_;
};

}  // namespace mbse
