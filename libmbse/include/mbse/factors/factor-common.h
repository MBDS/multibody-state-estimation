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

#include <gtsam/base/Vector.h>
#include <gtsam/base/VectorSpace.h>

namespace mbse
{
/** Type for system internal states q_{k}, dq_{k}, ddq_{k} */
class state_t : public gtsam::Vector
{
   public:
	state_t() = default;
	state_t(double val) : gtsam::Vector(1) { (*this)[0] = val; }
	state_t(const gtsam::Vector& v) : gtsam::Vector(v) {}

	// --- Interface expected by gtsam::VectorSpace ---
	/** print with optional string */
	void print(const std::string& prefix = "") const;
	/** equals with an tolerance */
	bool equals(const state_t& p, double tol = 1e-9) const
	{
		return gtsam::traits<gtsam::Vector>::Equals(*this, p, tol);
	}
	enum
	{
		dimension = -1
	};
	int dim() const { return gtsam::Vector::size(); }
	/// identity for group operation
	inline static state_t identity() { return state_t(0.0); }
	/// return as Eigen::Vector
	const gtsam::Vector& vector() const { return *this; }

	template <typename vector_like_t>
	state_t operator+(const vector_like_t& b) const
	{
		state_t a = *this;
		a += b;
		return a;
	}
	state_t operator-(const state_t& b) const
	{
		state_t a = *this;
		a -= b;
		return a;
	}
	state_t operator-() const
	{
		state_t a = *this;
		a *= -1;
		return a;
	}
	// --- End of interface expected by gtsam::VectorSpace ---
};

}  // namespace mbse

namespace gtsam
{
// traits for: state_t
template <>
struct traits<mbse::state_t> : public internal::VectorSpace<mbse::state_t>
{
};
template <>
struct traits<const mbse::state_t> : public internal::VectorSpace<mbse::state_t>
{
};

}  // namespace gtsam
