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

#include <mbse/factors/factor-common.h>
#include <mbse/dynamics/dynamic-simulators.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

namespace mbse
{
/** Factor for multibody forward (direct) dynamics in dependent coordinates,
 * using a general dynamics solver.
 *
 * This implements: \f$ error = \ddot{q}_k -
 * \ddot{q}(q_k,\dot{q}_k,M,...) = 0\f$
 *
 * Unknowns: \f$ q_k, \dot{q}_k, \ddot{q}_k \f$
 *
 * Fixed data: multibody model (inertias, masses, etc.), external forces.
 */
class FactorDynamics
	: public gtsam::NoiseModelFactor3<
		  state_t /* q_k */, state_t /* dq_k */, state_t /* ddq_k */>
{
   private:
	using This = FactorDynamics;
	using Base = gtsam::NoiseModelFactor3<state_t, state_t, state_t>;

	/** This factor's own private, exclusively-owned clone of the solver
	 * (and, transitively, the AssembledRigidModel) passed in at construction.
	 * Owning a private clone -instead of sharing the caller's instance-
	 * is what makes it safe to evaluate different factors concurrently, e.g.
	 * from different threads as GTSAM does when built with TBB. */
	CDynamicSimulatorBase::Ptr dynamic_solver_owned_;
	CDynamicSimulatorBase* dynamic_solver_ = nullptr;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorDynamics() = default;

	/** Constructor. Note that `dynamic_solver` is only used as a template to
	 * clone from: this factor does NOT keep a reference to it, nor mutates it.
	 */
	FactorDynamics(
		CDynamicSimulatorBase* dynamic_solver,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dq_k, gtsam::Key key_ddq_k)
		: Base(noiseModel, key_q_k, key_dq_k, key_ddq_k),
		  dynamic_solver_owned_(dynamic_solver->clone()),
		  dynamic_solver_(dynamic_solver_owned_.get())
	{
	}

	FactorDynamics(const FactorDynamics& o)
		: Base(o),
		  dynamic_solver_owned_(o.dynamic_solver_owned_->clone()),
		  dynamic_solver_(dynamic_solver_owned_.get())
	{
	}
	FactorDynamics& operator=(const FactorDynamics& o) = delete;

	virtual ~FactorDynamics() override;

	/// @return a deep copy of this factor
	virtual gtsam::NonlinearFactor::shared_ptr clone() const override;

	/** implement functions needed for Testable */

	/** print */
	virtual void print(
		const std::string& s, const gtsam::KeyFormatter& keyFormatter =
								  gtsam::DefaultKeyFormatter) const override;

	/** equals */
	virtual bool equals(
		const gtsam::NonlinearFactor& expected,
		double tol = 1e-9) const override;

	/** implement functions needed to derive from Factor */

	/** vector of errors */
	gtsam::Vector evaluateError(
		const state_t& q_k, const state_t& dq_k, const state_t& ddq_k,
		MBSE_OptionalMatrixType H1 = OptionalNone,
		MBSE_OptionalMatrixType H2 = OptionalNone,
		MBSE_OptionalMatrixType H3 = OptionalNone) const override;

	/** number of variables attached to this factor */
	std::size_t size() const { return 3; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
		ar& boost::serialization::make_nvp(
			"FactorDynamics", boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
