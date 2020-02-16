
#pragma once

#include <sparsembs/factor-common.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <sparsembs/CAssembledRigidModel.h>

namespace sparsembs
{
/** Factor for contraints
 * Here the constraints equation Phi(q)=0 is implmented
 */

// Create derived class "ConstraintsFactor" from superclass "NoiseModelFactor3"
class FactorConstraints : public gtsam::NoiseModelFactor1<state_t>
{
   private:
	using This = FactorConstraints;
	using Base = gtsam::NoiseModelFactor1<state_t>;

	// Class parameters (pointer to type "CConstraintBase")
	CAssembledRigidModel::Ptr m_arm;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraints() = default;

	/** Construcotr */
	FactorConstraints(
		const CAssembledRigidModel::Ptr& arm,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k)
		: Base(noiseModel, key_q_k), m_arm(arm)
	{
	}

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
		boost::optional<gtsam::Matrix&> H1 = boost::none) const override;

	/** numberof variable attached to this factor */
	std::size_t size() const { return 1; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
		ar& boost::serialization::make_nvp(
			"FactorConstraints",
			boost::serialization::base_object<Base>(*this));
	}
};

}  // namespace sparsembs