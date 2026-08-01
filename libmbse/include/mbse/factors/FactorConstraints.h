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
/** Factor for generalized coordinate contraints Phi(q)=0.
 */
class FactorConstraints : public gtsam::NoiseModelFactor1<state_t>
{
   private:
	using This = FactorConstraints;
	using Base = gtsam::NoiseModelFactor1<state_t>;

	/** This factor's own private, exclusively-owned clone of the model passed
	 * in at construction. Owning a private clone -instead of sharing the
	 * caller's instance- is what makes it safe to evaluate different factors
	 * concurrently, e.g. from different threads as GTSAM does when built with
	 * TBB. */
	AssembledRigidModel::Ptr arm_;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraints() = default;

	/** Constructor. Note that `arm` is only used as a template to clone from:
	 * this factor does NOT keep a reference to it, nor mutates it. */
	FactorConstraints(
		const AssembledRigidModel::Ptr& arm,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k)
		: Base(noiseModel, key_q_k), arm_(arm->clone())
	{
	}

	FactorConstraints(const FactorConstraints& o)
		: Base(o), arm_(o.arm_->clone())
	{
	}
	FactorConstraints& operator=(const FactorConstraints& o) = delete;

	virtual ~FactorConstraints() override;

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
		const state_t& q_k,
		MBSE_OptionalMatrixType H1 = OptionalNone) const override;

	/** numberof variable attached to this factor */
	std::size_t size() const { return 1; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
#ifdef GTSAM_ENABLE_BOOST_SERIALIZATION
		ar& boost::serialization::make_nvp(
			"FactorConstraints",
			boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
