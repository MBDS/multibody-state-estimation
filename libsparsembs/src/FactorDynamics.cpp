#include <sparsembs/FactorDynamics.h>
#include <sparsembs/CAssembledRigidModel.h>

using namespace sparsembs;

FactorDynamics::~FactorDynamics() = default;

gtsam::NonlinearFactor::shared_ptr FactorDynamics::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorDynamics::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorDynamics(" << keyFormatter(this->key1()) << ","
			  << keyFormatter(this->key2()) << "," << keyFormatter(this->key3())
			  << ")\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
}

bool FactorDynamics::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   m_dynamic_solver->get_model() == e->m_dynamic_solver->get_model();
}

gtsam::Vector FactorDynamics::evaluateError(
	const state_t& q_k, const state_t& dq_k, const state_t& ddq_k,
	boost::optional<gtsam::Matrix&> H1, boost::optional<gtsam::Matrix&> H2,
	boost::optional<gtsam::Matrix&> H3) const
{
	const auto n = q_k.size();
	if (dq_k.size() != n || ddq_k.size() != n)
		throw std::runtime_error("Inconsistent vector lengths!");
	if (n < 1) throw std::runtime_error("Empty state vector!");

	// Set q & dq in the multibody model:
	CAssembledRigidModel& arm = *m_dynamic_solver->get_model_non_const();
	arm.m_q = q_k.vector();
	arm.m_dotq = dq_k.vector();

	// Predict accelerations:
	Eigen::VectorXd qpp_predicted;
    const double t = 0;  // wallclock time (useless?)

	m_dynamic_solver->solve_ddotq(t, qpp_predicted);

	// Evaluate error:
	gtsam::Vector err = qpp_predicted - ddq_k.vector();

	// d err / d q_k
	if (H1)
	{
		auto& Hv = H1.value();
		Hv.setZero(n, n);
	}
	// d err / d dq_k
	if (H2)
	{
		auto& Hv = H2.value();
		Hv.setZero(n, n);
	}
	// d err / d ddq_k
	if (H3)
	{
		auto& Hv = H3.value();
		Hv = -Eigen::MatrixXd::Identity(n, n);
	}
	return err;
}
