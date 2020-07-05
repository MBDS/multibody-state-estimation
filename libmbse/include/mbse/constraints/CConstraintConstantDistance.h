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
/** Constraint: constant distance between two points */
class CConstraintConstantDistance : public CConstraintBase,
									public CConstraintCommon<2>
{
   public:
	using me_t = CConstraintConstantDistance;

	double length;

	CConstraintConstantDistance(
		const size_t _point_index0, const size_t _point_index1,
		const double _length)
		: CConstraintCommon({_point_index0, _point_index1}), length(_length)
	{
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }
};

}  // namespace mbse
