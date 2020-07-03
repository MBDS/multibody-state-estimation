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

#include <mbse/CConstraintBase.h>
#include <mbse/CConstraintCommon.h>

namespace mbse
{
/** Constraint (for relative coordinates): relative angle in "point 1" between
 * two rods pt1-pt2 and pt1-pt3 (CCW=positive):
 * \code
 *        pt2
 *       /
 *      /
 *     / ) angle
 * pt0 ----------- pt1
 * \endcode
 */
class CConstraintRelativeAngle : public CConstraintBase,
								 public CConstraintCommon<3>
{
   public:
	using me_t = CConstraintRelativeAngle;

	size_t angleIndexInQ = static_cast<size_t>(-1);

	CConstraintRelativeAngle(
		const size_t _point_index0, const size_t _point_index1,
		const size_t _point_index2, const size_t _angleIndexInQ)
		: CConstraintCommon({_point_index0, _point_index1, _point_index2}),
		  angleIndexInQ(_angleIndexInQ)
	{
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }

   protected:
};

}  // namespace mbse
