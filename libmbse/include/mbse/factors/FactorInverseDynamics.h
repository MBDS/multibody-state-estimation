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
/** Factor for multibody inverse dynamics, using a general dynamics
 * solver, and given a list of actuated degrees of freedom (a subset of all
 * generized coordinates "q").
 *
 * This implements: \f$ error = \ddot{q}_k -
 * \ddot{q}(q_k,\dot{q}_k,Q_k) = 0\f$
 *
 * Unknowns: \f$ Q_k \f$
 * Assumed to be fixed: \f$ q_k, \dot{q}_k, \ddot{q}_k \f$
 *
 * Note that this factor only provides the Jacobian wrt the external
 * forces Q_k, disregarding the Jacobians wrt the other variables.
 * This is done in order to leave at the factor graph designer's criterion
 * how velocities and positions are integrated over time, from which the implied
 * forces are obtained via this factor without modifying trajectories.
 *
 * Fixed data: multibody model (inertias, masses, etc.)
 */
class FactorInverseDynamics : public gtsam::NoiseModelFactor1<state_t /* Q_k */>
{
   private:
	using This = FactorInverseDynamics;
	using Base = gtsam::NoiseModelFactor1<state_t>;

	CDynamicSimulatorBase* dynamic_solver_ = nullptr;

	gtsam::Key key_q_k_, key_dq_k_, key_ddq_k_;
	const gtsam::Values* valuesFor_q_dq_ = nullptr;
	std::vector<size_t> actuatedDegreesOfFreedomInQ_;

	mutable state_t cached_q_, cached_dq_, cached_ddq_, cached_Q_;
	mutable gtsam::Matrix cached_d_e_Q_;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorInverseDynamics() = default;

	/** Constructor */
	FactorInverseDynamics(
		CDynamicSimulatorBase* dynamic_solver,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dq_k, gtsam::Key key_ddq_k, gtsam::Key key_Q_k,
		const gtsam::Values& valuesFor_q_dq,
		const std::vector<size_t>& actuatedDegreesOfFreedomInQ)
		: Base(noiseModel, key_Q_k),
		  dynamic_solver_(dynamic_solver),
		  key_q_k_(key_q_k),
		  key_dq_k_(key_dq_k),
		  key_ddq_k_(key_ddq_k),
		  valuesFor_q_dq_(&valuesFor_q_dq),
		  actuatedDegreesOfFreedomInQ_(actuatedDegreesOfFreedomInQ)
	{
	}

	virtual ~FactorInverseDynamics() override;

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
		const state_t& Q_k,
		gtsam::OptionalMatrixType d_e_Q = OptionalNone) const override;

	/** number of variables attached to this factor */
	std::size_t size() const { return 2; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
		ar& boost::serialization::make_nvp(
			"FactorInverseDynamics",
			boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
