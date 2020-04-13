/*+-------------------------------------------------------------------------+
  |    XXX
  |                                                                         |
  | Copyright (C) 2019-2020 University of Almeria                           |
  | See README for list of authors and papers                               |
  | Distributed under GNU General Public License version 3                  |
  |   See <http://www.gnu.org/licenses/>                                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/factor-common.h>
#include <mbse/dynamic-simulators.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

namespace mbse
{
/** Factor for multibody dynamics, using a general dynamics solver.
 *
 * This implements: \f$ \text{error} = \ddot{q}_k -
 * \ddot{q}(q_k,\dot{q}_k,M,...) = 0\f$
 *
 * Unknowns: \ddot{q}_k, \dot{q}_k, q_k
 *
 * Fixed data:
 */
class FactorDynamics
	: public gtsam::NoiseModelFactor3<state_t, state_t, state_t>
{
   private:
	using This = FactorDynamics;
	using Base = gtsam::NoiseModelFactor3<state_t, state_t, state_t>;

	CDynamicSimulatorBase* m_dynamic_solver = nullptr;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorDynamics() = default;

	/** Constructor */
	FactorDynamics(
		CDynamicSimulatorBase* dynamic_solver,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dq_k, gtsam::Key key_ddq_k)
		: Base(noiseModel, key_q_k, key_dq_k, key_ddq_k),
		  m_dynamic_solver(dynamic_solver)
	{
	}

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
		boost::optional<gtsam::Matrix&> H1 = boost::none,
		boost::optional<gtsam::Matrix&> H2 = boost::none,
		boost::optional<gtsam::Matrix&> H3 = boost::none) const override;

	/** number of variables attached to this factor */
	std::size_t size() const { return 3; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
		ar& boost::serialization::make_nvp(
			"FactorDynamics", boost::serialization::base_object<Base>(*this));
	}
};

}  // namespace mbse
