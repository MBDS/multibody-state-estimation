/*+-------------------------------------------------------------------------+
  |    FactorGraph Control (sparsembs)  C++ library                         |
  |                                                                         |
  | Copyright (C) 2019-2020 University of Almeria                           |
  | See README for list of authors and papers                               |
  | Distributed under GNU General Public License version 3                  |
  |   See <http://www.gnu.org/licenses/>                                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/FactorTrapInt.h>

using namespace mbse;

FactorTrapInt::~FactorTrapInt() = default;

gtsam::NonlinearFactor::shared_ptr FactorTrapInt::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

// Build function �print� defined in the header �FactorTrapInt.h�
void FactorTrapInt::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorTrapInt(" << keyFormatter(this->key1()) << ","
			  << keyFormatter(this->key2()) << "," << keyFormatter(this->key3())
			  << "," << keyFormatter(this->key4()) << ")\n";
	gtsam::traits<double>::Print(timestep_, "  timestep: ");
	this->noiseModel_->print("  noise model: ");
}

bool FactorTrapInt::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   gtsam::traits<double>::Equals(this->timestep_, e->timestep_, tol);
}

// Build function �evaluateError� defined in the header �FactorTrapInt.h�

gtsam::Vector FactorTrapInt::evaluateError(
	const state_t& x_k, const state_t& x_kp1, const state_t& v_k,
	const state_t& v_kp1, boost::optional<gtsam::Matrix&> H1,
	boost::optional<gtsam::Matrix&> H2, boost::optional<gtsam::Matrix&> H3,
	boost::optional<gtsam::Matrix&> H4) const
{
	const auto n = x_k.size();
	if (x_kp1.size() != n || v_k.size() != n)
		throw std::runtime_error("Inconsistent vector lengths!");

	gtsam::Vector err = x_kp1.vector() - x_k.vector() -
						0.5 * timestep_ * v_k.vector() -
						0.5 * timestep_ * v_kp1.vector();

	// Jacobian of �err� respect to[x_k x_kp1 v_k v_kp1]
	if (H1)
	{
		auto& H1v = H1.value();
		H1v = -Eigen::MatrixXd::Identity(n, n);
	}
	if (H2)
	{
		auto& H2v = H2.value();
		H2v = Eigen::MatrixXd::Identity(n, n);
	}
	if (H3)
	{
		auto& H3v = H3.value();
		H3v = -0.5 * timestep_ * Eigen::MatrixXd::Identity(n, n);
	}
	if (H4)
	{
		auto& H4v = H4.value();
		H4v = -0.5 * timestep_ * Eigen::MatrixXd::Identity(n, n);
	}

	return err;
}
