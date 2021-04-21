/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factors/FactorConstraintsVelIndep.h>
#include <mbse/AssembledRigidModel.h>
#include <mrpt/core/exceptions.h>
#include <mbse/mbse-utils.h>

using namespace mbse;

FactorConstraintsVelIndep::FactorConstraintsVelIndep(
	const AssembledRigidModel::Ptr& arm,
	const std::vector<size_t>& indCoordsIndices,
	const gtsam::SharedNoiseModel& noiseModel, gtsam::Key key_q_k,
	gtsam::Key key_dotq_k, gtsam::Key key_dotz_k)
	: Base(noiseModel, key_q_k, key_dotq_k, key_dotz_k),
	  arm_(arm),
	  indCoordsIndices_(indCoordsIndices),
	  matrix_Iidx_(selector_matrix(indCoordsIndices, arm->q_.size()))
{
}

FactorConstraintsVelIndep::~FactorConstraintsVelIndep() = default;

gtsam::NonlinearFactor::shared_ptr FactorConstraintsVelIndep::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorConstraintsVelIndep::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "mbse::FactorConstraintsVelIndep("
			  << keyFormatter(this->key1()) << "," << keyFormatter(this->key2())
			  << "," << keyFormatter(this->key3()) << ")\n";
	noiseModel_->print("  noise model: ");
}

bool FactorConstraintsVelIndep::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol);
}

gtsam::Vector FactorConstraintsVelIndep::evaluateError(
	const state_t& q_k, const state_t& dotq_k, const state_t& dotz_k,
	boost::optional<gtsam::Matrix&> de_dq,
	boost::optional<gtsam::Matrix&> de_dqp,
	boost::optional<gtsam::Matrix&> de_dzp) const
{
	MRPT_START

	const auto n = q_k.size();
	if (n < 1) throw std::runtime_error("Empty state vector q_k!");
	const auto d = dotz_k.size();
	if (d < 1) throw std::runtime_error("Empty state vector dotz_k!");

	ASSERT_EQUAL_(dotq_k.size(), q_k.size());
	ASSERT_(q_k.size() > 0);

	// Set q in the multibody model:
	arm_->q_ = q_k;
	arm_->dotq_ = dotq_k;

	// Update Jacobian and Hessian tensor:
	arm_->update_numeric_Phi_and_Jacobians();

	const auto m = arm_->Phi_.rows();
	if (m < 1) throw std::runtime_error("Empty Phi() vector!");

	// Evaluate error:
	const Eigen::MatrixXd Phi_q = arm_->Phi_q_.asDense();
	const Eigen::MatrixXd dPhiqdq_dq = Phi_q;

	gtsam::Vector err = gtsam::Vector::Zero(m + d);
	err.head(m) = Phi_q * dotq_k;
	err.tail(d) = mbse::subset(dotq_k, indCoordsIndices_) - dotz_k;

	// Get the Jacobians required for optimization:
	// (Section 6.8 of the paper)
	// d err / d q_k
	if (de_dq)
	{
		auto& Hv = de_dq.value();
		Hv.resize(m + d, n);
		// Phi_qq*dq = \dot{Phi_q}
		Hv.block(0, 0, m, n) = arm_->dotPhi_q_.asDense();
		Hv.block(m, 0, d, n).setZero();
	}

	if (de_dqp)
	{
		auto& Hv = de_dqp.value();
		Hv.resize(m + d, n);
		Hv.block(0, 0, m, n) = arm_->Phi_q_.asDense();
		Hv.block(m, 0, d, n) = matrix_Iidx_;
	}

	if (de_dzp)
	{
		auto& Hv = de_dzp.value();
		Hv.setZero(m + d, d);
		Hv.block(m, 0, d, d) = -gtsam::Matrix::Identity(d, d);
	}

	return err;

	MRPT_END
}
