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

#include <gtsam/geometry/Point3.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

namespace mbse
{
/** Factor for numerical integration using the Euler methods.
 *
 * This implements: \f$x_{k+1} = x_{k} + dt * v_{k}\f$
 *
 * Unknowns: x_{k},x_{k+1}, v_{k}
 * Fixed data: dt
 */
class FactorEulerInt
	: public gtsam::NoiseModelFactor3<state_t, state_t, state_t>
{
   private:
	using This = FactorEulerInt;
	using Base = gtsam::NoiseModelFactor3<state_t, state_t, state_t>;

	/** Numerical integration timestep */
	double timestep_ = 0;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorEulerInt() = default;

	/** Constructor */
	FactorEulerInt(
		const double timestep, const gtsam::SharedNoiseModel& noiseModel,
		gtsam::Key key_x_k, gtsam::Key key_x_kp1, gtsam::Key key_v_k)
		: Base(noiseModel, key_x_k, key_x_kp1, key_v_k), timestep_(timestep)
	{
	}

	virtual ~FactorEulerInt() override;

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
		const state_t& x_k, const state_t& x_kp1, const state_t& v_k,
		gtsam::OptionalMatrixType H1 = OptionalNone,
		gtsam::OptionalMatrixType H2 = OptionalNone,
		gtsam::OptionalMatrixType H3 = OptionalNone) const override;

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
			"FactorEulerInt", boost::serialization::base_object<Base>(*this));
		ar& BOOST_SERIALIZATION_NVP(timestep_);
#endif
	}
};

}  // namespace mbse
