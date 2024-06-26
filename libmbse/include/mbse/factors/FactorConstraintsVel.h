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
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <mbse/AssembledRigidModel.h>

namespace mbse
{
/** Factor for contraints velocity
 * Here the constraints equation Phiq(q)*dotq=0 is implemented; for fixed
 * constraints, Phit(q)=0.
 */
class FactorConstraintsVel : public gtsam::NoiseModelFactor2<state_t, state_t>
{
   private:
	using This = FactorConstraintsVel;
	using Base = gtsam::NoiseModelFactor2<state_t, state_t>;

	// Class parameters (pointer to type "ConstraintBase")
	AssembledRigidModel::Ptr arm_;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraintsVel() = default;

	/** Construcotr */
	FactorConstraintsVel(
		const AssembledRigidModel::Ptr& arm,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dotq_k)
		: Base(noiseModel, key_q_k, key_dotq_k), arm_(arm)
	{
	}

	virtual ~FactorConstraintsVel() override;

	// @return a deep copy of this factor
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
		const state_t& q_k, const state_t& dotq_k,
		gtsam::OptionalMatrixType H1 = OptionalNone,
		gtsam::OptionalMatrixType H2 = OptionalNone) const override;

	/** numberof variable attached to this factor */
	std::size_t size() const { return 2; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
		ar& boost::serialization::make_nvp(
			"FactorConstraintsVel",
			boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
