/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#pragma once

#include <mbse/factors/factor-common.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <mbse/CAssembledRigidModel.h>

namespace mbse
{
/** Factor for contraints velocity (independent coordinates).
 * Here the constraints equation Phiq(q)*dotq=0 is implemented; for fixed
 * constraints, Phit(q)=0.
 */
class FactorConstraintsVelIndep
	: public gtsam::NoiseModelFactor3<
		  state_t /*Q*/, state_t /*dotQ*/, state_t /*dotZ*/>
{
   private:
	using This = FactorConstraintsVelIndep;
	using Base = gtsam::NoiseModelFactor3<
		state_t /*Q*/, state_t /*dotQ*/, state_t /*dotZ*/>;

	// Class parameters (pointer to type "CConstraintBase")
	CAssembledRigidModel::Ptr arm_;
	std::vector<size_t> indCoordsIndices_;
	gtsam::Matrix matrix_Iidx_;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = boost::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraintsVelIndep() = default;

	/** Construcotr */
	FactorConstraintsVelIndep(
		const CAssembledRigidModel::Ptr& arm,
		const std::vector<size_t>& indCoordsIndices,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
		gtsam::Key key_dotq_k, gtsam::Key key_z_k);

	virtual ~FactorConstraintsVelIndep() override;

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
		const state_t& q_k, const state_t& dotq_k, const state_t& dotz_k,
		boost::optional<gtsam::Matrix&> de_dq = boost::none,
		boost::optional<gtsam::Matrix&> de_dqp = boost::none,
		boost::optional<gtsam::Matrix&> de_dzp = boost::none) const override;

	/** numberof variable attached to this factor */
	std::size_t size() const { return 3; }

   private:
	/** Serialization function */
	friend class boost::serialization::access;
	template <class ARCHIVE>
	void serialize(ARCHIVE& ar, const unsigned int /*version*/)
	{
		ar& boost::serialization::make_nvp(
			"FactorConstraintsVelIndep",
			boost::serialization::base_object<Base>(*this));
	}
};

}  // namespace mbse
