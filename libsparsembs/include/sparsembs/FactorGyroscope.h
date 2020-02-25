/*+-------------------------------------------------------------------------+
  |    XXX
  |                                                                         |
  | Copyright (C) 2019-2020 University of Almeria                           |
  | See README for list of authors and papers                               |
  | Distributed under GNU General Public License version 3                  |
  |   See <http://www.gnu.org/licenses/>                                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <sparsembs/factor-common.h>
#include <sparsembs/dynamic-simulators.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <sparsembs/virtual-sensors.h>

namespace sparsembs
{
/** Factor for gyroscope reading.
 *
 * Fixed data:
 */
class FactorGyroscope : public gtsam::NoiseModelFactor2<state_t, state_t>
{
   private:
    using This = FactorGyroscope;
    using Base = gtsam::NoiseModelFactor2<state_t, state_t>;

    CVirtualSensor_Gyro* m_gyro = nullptr;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorGyroscope() = default;

	/** Constructor */
	FactorGyroscope(
		CVirtualSensor_Gyro* gyro_solver,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dq_k)
		: Base(noiseModel, key_q_k, key_dq_k), m_gyro(gyro_solver)
	{
	}

	virtual ~FactorGyroscope() override;

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
		boost::optional<gtsam::Matrix&> H1 = boost::none,
		boost::optional<gtsam::Matrix&> H2 = boost::none) const override;

	/** number of variables attached to this factor */
	std::size_t size() const { return 2; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
		ar& boost::serialization::make_nvp(
			"FactorGyroscope", boost::serialization::base_object<Base>(*this));
	}
};

}  // namespace sparsembs
