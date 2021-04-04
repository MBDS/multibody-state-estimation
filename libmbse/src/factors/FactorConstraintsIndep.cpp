/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factors/FactorConstraintsIndep.h>
#include <mbse/CAssembledRigidModel.h>
#include <mbse/mbse-utils.h>

using namespace mbse;

FactorConstraintsIndep::FactorConstraintsIndep(
	const CAssembledRigidModel::Ptr& arm,
	const std::vector<size_t>& indCoordsIndices,
	const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_z_k,
	gtsam::Key key_q_k)
	: Base(noiseModel, key_z_k, key_q_k),
	  arm_(arm),
	  indCoordsIndices_(indCoordsIndices),
	  matrix_Iidx_(selector_matrix(indCoordsIndices, arm->q_.size()))
{
}

FactorConstraintsIndep::~FactorConstraintsIndep() = default;

gtsam::NonlinearFactor::shared_ptr FactorConstraintsIndep::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorConstraintsIndep::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorConstraintsIndep("
			  << keyFormatter(this->key1()) << "," << keyFormatter(this->key2())
			  << ")\n";
	noiseModel_->print("  noise model: ");
}

bool FactorConstraintsIndep::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol);
}

gtsam::Vector FactorConstraintsIndep::evaluateError(
	const state_t& z_k, const state_t& q_k,
	boost::optional<gtsam::Matrix&> de_dz,
	boost::optional<gtsam::Matrix&> de_dq) const
{
	MRPT_START

	const auto n = q_k.size();
	if (n < 1) throw std::runtime_error("Empty state vector q_k!");
	const auto d = z_k.size();
	if (d < 1) throw std::runtime_error("Empty state vector z_k!");

	// Set q in the multibody model:
	arm_->q_ = q_k;

	// Update Jacobians:
	arm_->update_numeric_Phi_and_Jacobians();

	const auto m = arm_->Phi_.rows();
	if (m < 1) throw std::runtime_error("Empty Phi() vector!");

	// Evaluate error:
	gtsam::Vector err = gtsam::Vector::Zero(m + d);
	err.head(m) = arm_->Phi_;
	err.tail(d) = mbse::subset(q_k, indCoordsIndices_) - z_k;

	// Get the Jacobians required for optimization:
	// d err / d z_k
	if (de_dz)
	{
		auto& Hv = de_dz.value();
		Hv.setZero(m + d, d);
		Hv.block(m, 0, d, d) = -gtsam::Matrix::Identity(d, d);
	}

	// d err / d q_k
	if (de_dq)
	{
		auto& Hv = de_dq.value();
		Hv.resize(m + d, n);
		Hv.block(0, 0, m, n) = arm_->Phi_q_.asDense();
		// Build "I_idx", as called in the paper (sect. 6.7)
		Hv.block(m, 0, d, n) = matrix_Iidx_;
	}

	return err;

	MRPT_END
}
