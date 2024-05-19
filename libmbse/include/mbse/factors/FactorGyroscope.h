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
#include <mbse/virtual-sensors.h>

namespace mbse
{
/** Factor for 2D gyroscope reading.
 *
 * Fixed data:
 */
class FactorGyroscope : public gtsam::NoiseModelFactor2<state_t, state_t>
{
   private:
	using This = FactorGyroscope;
	using Base = gtsam::NoiseModelFactor2<state_t, state_t>;

	AssembledRigidModel* arm_ = nullptr;
	size_t body_idx_ = 0;
	double reading_ = 0;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorGyroscope() = default;
	virtual ~FactorGyroscope() override = default;

	/** Constructor. angvel_reading in rad/sec, positive CCW. */
	FactorGyroscope(
		AssembledRigidModel& arm, const size_t body_idx,
		const double angvel_reading, const gtsam::SharedNoiseModel& noiseModel,
		gtsam::Key key_q_k, gtsam::Key key_dq_k)
		: Base(noiseModel, key_q_k, key_dq_k),
		  arm_(&arm),
		  body_idx_(body_idx),
		  reading_(angvel_reading)
	{
	}

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
		const state_t& q_k, const state_t& dq_k,
		boost::optional<gtsam::Matrix&>  H1 = boost::none,
		boost::optional<gtsam::Matrix&>  H2 = boost::none) const override;

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
			"FactorGyroscope", boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
