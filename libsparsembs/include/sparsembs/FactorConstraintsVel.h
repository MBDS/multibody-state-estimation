
#pragma once

#include <sparsembs/factor-common.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <sparsembs/CAssembledRigidModel.h>

namespace sparsembs
{
/** Factor for contraints velocity
 * Here the constraints equation Phiq(q)*dotq=0 is implmented; for fixed
 * constraints, Phit(q)=0.
 */

// Create derived class "ConstraintsFactor" from superclass "NoiseModelFactor2"
class FactorConstraintsVel : public gtsam::NoiseModelFactor2<state_t, state_t>
{
   private:
    using This = FactorConstraintsVel;
    using Base = gtsam::NoiseModelFactor2<state_t, state_t>;

	// Class parameters (pointer to type "CConstraintBase")
	CAssembledRigidModel::Ptr m_arm;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraintsVel() = default;

	/** Construcotr */
	FactorConstraintsVel(
		const CAssembledRigidModel::Ptr& arm,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dotq_k)
		: Base(noiseModel, key_q_k, key_dotq_k), m_arm(arm)
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
		boost::optional<gtsam::Matrix&> H1 = boost::none,
		boost::optional<Eigen::Tensor<double, 3>&> H2 =
			boost::none) const override;
	//==========================================================================
	//  H2 is a tensor, but it seems to not exist in Eigen! Maybe Eigen version
	//  is not enough updated or I have to define tensorse before its using!
	//==========================================================================

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
