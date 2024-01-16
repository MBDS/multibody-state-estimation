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
/** Constraint: constant distance between two points */
class ConstraintConstantDistance : public ConstraintBase,
								   public ConstraintCommon<2>
{
   public:
	using me_t = ConstraintConstantDistance;

	double length;

	ConstraintConstantDistance(
		const size_t _point_index0, const size_t _point_index1,
		const double _length)
		: ConstraintCommon({_point_index0, _point_index1}), length(_length)
	{
	}

	void buildSparseStructures(AssembledRigidModel& arm) const override;
	void update(AssembledRigidModel& arm) const override;
	void print(std::ostream& o) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }
};

}  // namespace mbse
