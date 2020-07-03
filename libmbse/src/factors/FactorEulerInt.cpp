/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/factors/FactorEulerInt.h>
#include <mrpt/core/exceptions.h>

using namespace mbse;

FactorEulerInt::~FactorEulerInt() = default;

gtsam::NonlinearFactor::shared_ptr FactorEulerInt::clone() const
{
	return boost::static_pointer_cast<gtsam::NonlinearFactor>(
		gtsam::NonlinearFactor::shared_ptr(new This(*this)));
}

void FactorEulerInt::print(
	const std::string& s, const gtsam::KeyFormatter& keyFormatter) const
{
	std::cout << s << "FactorEulerInt(" << keyFormatter(this->key1()) << ","
			  << keyFormatter(this->key2()) << "," << keyFormatter(this->key3())
			  << ")\n";
	gtsam::traits<double>::Print(timestep_, "  timestep: ");
	noiseModel_->print("  noise model: ");
}

bool FactorEulerInt::equals(
	const gtsam::NonlinearFactor& expected, double tol) const
{
	const This* e = dynamic_cast<const This*>(&expected);
	return e != nullptr && Base::equals(*e, tol) &&
		   gtsam::traits<double>::Equals(timestep_, e->timestep_, tol);
}

gtsam::Vector FactorEulerInt::evaluateError(
	const state_t& x_k, const state_t& x_kp1, const state_t& v_k,
	boost::optional<gtsam::Matrix&> H1, boost::optional<gtsam::Matrix&> H2,
	boost::optional<gtsam::Matrix&> H3) const
{
	const auto n = x_k.size();

	ASSERT_EQUAL_(x_kp1.size(), x_k.size());
	ASSERT_EQUAL_(v_k.size(), x_k.size());

	gtsam::Vector err =
		x_kp1.vector() - x_k.vector() - timestep_ * v_k.vector();

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
		H3v = -timestep_ * Eigen::MatrixXd::Identity(n, n);
	}

	return err;
}
