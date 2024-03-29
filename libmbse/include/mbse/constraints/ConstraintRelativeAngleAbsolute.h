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

#include <mbse/constraints/ConstraintBase.h>
#include <mbse/constraints/ConstraintCommon.h>

namespace mbse
{
/** Constraint (for relative coordinates): relative angle in "point 0" between
 * the absolute +X axis and the FIXED-LENGTH rod pt0-pt1 (CCW=positive):
 * \code
 *        pt1
 *       /
 *      /
 *     / ) angle
 * pt0 -------------> +X axis
 * \endcode
 */
class ConstraintRelativeAngleAbsolute
	: public ConstraintBase,
	  public ConstraintCommon<
		  2 /* Num Euclidean points*/, 1 /* Num relative coords */,
		  1 /* Num Jacobian rows*/>
{
   public:
	using me_t = ConstraintRelativeAngleAbsolute;

	ConstraintRelativeAngleAbsolute(
		const size_t _point_index0, const size_t _point_index1,
		const size_t _angleIndexInQ)
		: ConstraintCommon({_point_index0, _point_index1}, {_angleIndexInQ})
	{
	}

	void buildSparseStructures(AssembledRigidModel& arm) const override;
	void update(AssembledRigidModel& arm) const override;
	void print(std::ostream& o) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }

	void realizeOperatingPoint(const AssembledRigidModel& arm) const override;

   protected:
	/** Proxy for length between the two points */
	mutable double L_ = .0;
	mutable bool useCos_ = true;
};

}  // namespace mbse
