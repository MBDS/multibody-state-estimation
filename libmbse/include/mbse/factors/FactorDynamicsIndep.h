/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
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
/** Factor for multibody forward (direct) dynamics in independent coordinates,
 * using a general dynamics solver.
 *
 * This implements: \f$ error = \ddot{z}_k -
 * \ddot{z}(q_k,\dot{q}_k,M,...) = 0\f$
 *
 * Unknowns: \f$ z_k, \dot{z}_k, \ddot{z}_k \f$
 *
 * Fixed data: multibody model (inertias, masses, etc.), external forces.
 *
 * Independent coordinates indices are taken from
 * CDynamicSimulatorIndepBase::independent_coordinate_indices()
 *
 */
class FactorDynamicsIndep
	: public gtsam::NoiseModelFactor3<
		  state_t /* z_k */, state_t /* dz_k */, state_t /* ddz_k */>
{
   private:
	using This = FactorDynamicsIndep;
	using Base = gtsam::NoiseModelFactor3<state_t, state_t, state_t>;

	CDynamicSimulatorIndepBase* dynamic_solver_ = nullptr;
	gtsam::Key key_q_k_;
	const gtsam::Values* valuesForQk_ = nullptr;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorDynamicsIndep() = default;

	/** Constructor */
	FactorDynamicsIndep(
		CDynamicSimulatorIndepBase* dynamic_solver,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_z_k,
		gtsam::Key key_dz_k, gtsam::Key key_ddz_k, gtsam::Key key_q_k,
		const gtsam::Values& valuesForQk)
		: Base(noiseModel, key_z_k, key_dz_k, key_ddz_k),
		  dynamic_solver_(dynamic_solver),
		  key_q_k_(key_q_k),
		  valuesForQk_(&valuesForQk)
	{
		MRPT_START
		ASSERTMSG_(
			!dynamic_solver->can_choose_indep_coords_,
			"Passed `dynamic_solver` must be configured NOT to dynamically "
			"choose DOF z variables for use within this factor.");
		MRPT_END
	}

	virtual ~FactorDynamicsIndep() override;

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
		const state_t& z_k, const state_t& dz_k, const state_t& ddz_k,
		boost::optional<gtsam::Matrix&> de_dz = boost::none,
		boost::optional<gtsam::Matrix&> de_dzp = boost::none,
		boost::optional<gtsam::Matrix&> de_dzpp = boost::none) const override;

	/** number of variables attached to this factor */
	std::size_t size() const { return 3; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
		ar& boost::serialization::make_nvp(
			"FactorDynamicsIndep",
			boost::serialization::base_object<Base>(*this));
	}
};

}  // namespace mbse
