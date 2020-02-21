#include <sparsembs/FactorConstraintsVel.h>
#include <sparsembs/CAssembledRigidModel.h>

#include <unsupported/Eigen/CXX11/Tensor>

using namespace sparsembs;

FactorConstraintsVel::~FactorConstraintsVel() = default;

gtsam::NonlinearFactor::shared_ptr FactorConstraintsVel::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorConstraintsVel::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorConstrintsVel(" << keyFormatter(this->key1())
			  << "," << keyFormatter(this->key2()) << ")\n";
	this->noiseModel_->print("  noise model: ");
}

bool FactorConstraintsVel::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol);
}

gtsam::Vector FactorConstraintsVel::evaluateError(
	const state_t& q_k, const state_t& dotq_k,
	boost::optional<gtsam::Matrix&> H1,
	boost::optional<gtsam::Matrix&> H2) const
{
	const auto n = q_k.size();
	if (dotq_k.size() != n)
		throw std::runtime_error("Inconsistent vector lengths!");
	if (n < 1) throw std::runtime_error("Empty state vector!");

	// Set q in the multibody model:
	m_arm->m_q = q_k.vector();
	m_arm->m_dotq = dotq_k.vector();

	// Update Jacobian and Hessian tensor:
	m_arm->update_numeric_Phi_and_Jacobians();

	// Evaluate error:
	const Eigen::MatrixXd Phi_q = m_arm->getPhi_q_dense();

	gtsam::Vector err = Phi_q * dotq_k.vector();

	// Get the Jacobians required for optimization:
	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
		Hv = m_arm->getPhi_q_dense();
	}

	if (H2)
	{
		auto& Hv = H2.value();
		// Hv = m_arm->  // getPhi_qq_dense(); Phi_qq is the Hessian Tensor
	}

	return err;
}
