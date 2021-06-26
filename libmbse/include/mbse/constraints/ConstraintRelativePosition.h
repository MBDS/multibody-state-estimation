/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
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
/** Constraint: relative position of one point wrt a system of coordinates built
 * from other three points. The relative position is automatically calculated
 * from "q" upon the first call to buildSparseStructures().
 *
 * See: Alejo Avello's book, section 4.4.1
 */
class ConstraintRelativePosition
	: public ConstraintBase,
	  public ConstraintCommon<
		  4 /*num points*/, 0 /*rel coords*/, 2 /*num constraint rows*/>
{
   public:
	using me_t = ConstraintRelativePosition;

	double x_ = 0, y_ = 0;

	ConstraintRelativePosition(
		const size_t refPointIndex0, const size_t refPointIndex1,
		const size_t refPointIndex2, const size_t constrainedPointIndex)
		: ConstraintCommon(
			  {refPointIndex0, refPointIndex1, refPointIndex2,
			   constrainedPointIndex})
	{
	}

	void buildSparseStructures(AssembledRigidModel& arm) const override;
	void update(AssembledRigidModel& arm) const override;
	void print(std::ostream& o) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }
};

}  // namespace mbse
