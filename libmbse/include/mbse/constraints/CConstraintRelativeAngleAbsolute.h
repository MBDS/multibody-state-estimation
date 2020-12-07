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

#include <mbse/constraints/CConstraintBase.h>
#include <mbse/constraints/CConstraintCommon.h>

namespace mbse
{
/** Constraint (for relative coordinates): relative angle in "point 0" between
 * the absolute +X axis and the rod pt0-pt1 (CCW=positive):
 * \code
 *        pt1
 *       /
 *      /
 *     / ) angle
 * pt0 -------------> +X axis
 * \endcode
 */
class CConstraintRelativeAngleAbsolute
	: public CConstraintBase,
	  public CConstraintCommon<
		  2 /* Num Euclidean points*/, 1 /* Num relative coords */,
		  1 /* Num Jacobian rows*/>
{
   public:
	using me_t = CConstraintRelativeAngleAbsolute;

	CConstraintRelativeAngleAbsolute(
		const size_t _point_index0, const size_t _point_index1,
		const size_t _angleIndexInQ)
		: CConstraintCommon({_point_index0, _point_index1}, {_angleIndexInQ})
	{
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }

   protected:
	/** Proxy for length between the two points */
	mutable double L_ = .0, Lsqr_ = .0;
};

}  // namespace mbse
