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
/** Factor for generalized *independent* coordinate contraints Phi(q)=0.
 */
class FactorConstraintsIndep
	: public gtsam::NoiseModelFactor2<state_t /*z_t*/, state_t /*q_t*/>
{
   private:
	using This = FactorConstraintsIndep;
	using Base = gtsam::NoiseModelFactor2<state_t /*z_t*/, state_t /*q_t*/>;

	// Class parameters (pointer to type "ConstraintBase")
	AssembledRigidModel::Ptr arm_;
	std::vector<size_t> indCoordsIndices_;
	gtsam::Matrix matrix_Iidx_;

   public:
	// shorthand for a smart pointer to a factor
	using shared_ptr = std::shared_ptr<This>;

	/** default constructor - only use for serialization */
	FactorConstraintsIndep() = default;

	/** Construcotr */
	FactorConstraintsIndep(
		const AssembledRigidModel::Ptr& arm,
		const std::vector<size_t>& indCoordsIndices,
		const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_z_k,
		gtsam::Key key_q_k);

	virtual ~FactorConstraintsIndep() override;

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
		const state_t& z_k, const state_t& q_k,
		gtsam::OptionalMatrixType de_dz = OptionalNone,
		gtsam::OptionalMatrixType de_dq = OptionalNone) const override;

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
			"FactorConstraintsIndep",
			boost::serialization::base_object<Base>(*this));
#endif
	}
};

}  // namespace mbse
