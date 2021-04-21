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

#include <mbse/constraints/ConstraintBase.h>
#include <mbse/constraints/ConstraintCommon.h>

namespace mbse
{
/** Constraint: forces a point to lie exactly on the line defined by two other
 * reference points */
class ConstraintMobileSlider : public ConstraintBase,
								public ConstraintCommon<3>
{
   public:
	using me_t = ConstraintMobileSlider;

	// Point indices: 0=mobile at the slider, 1 & 2: define the line on which 0
	// slides.

	ConstraintMobileSlider(
		const size_t _point_index, const size_t _ref_pt0, const size_t _ref_pt1)
		: ConstraintCommon({_point_index, _ref_pt0, _ref_pt1})
	{
	}

	void buildSparseStructures(AssembledRigidModel& arm) const override;
	void update(AssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }
};

}  // namespace mbse
