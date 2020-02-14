#include <sparsembs/FactorConstraints.h>
#include <sparsembs/CAssembledRigidModel.h>

using namespace sparsembs;

FactorConstraints::~FactorConstraints() = default;

gtsam::NonlinearFactor::shared_ptr FactorConstraints::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorConstraints::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorConstrints(" << keyFormatter(this->key()) << ")\n";
	// gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
}

/*bool FactorConstraints::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   m_constraints->get_mode () == e->m_constraints->get_model();
}*/

gtsam::Vector FactorConstraints::evaluateError(
	const state_t& q_k, boost::optional<gtsam::Matrix&> H1) const
{
	const auto n = q_k.size();
	if (n < 1) throw std::runtime_error("Empty state vector!");

	// Set q in the multibody model:
	CAssembledRigidModel& arm = *m_constraints->buildSparseStructures();
	arm.m_q = q_k.vector();
	CAssembledRigidModel& arm = *m_constraints->update();
	arm.m_Phi_q = H1;

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
