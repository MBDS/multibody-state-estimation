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
/** Constraint: forces a point to lie exactly on the line defined by two other
 * reference points */
class CConstraintMobileSlider : public CConstraintBase,
								public CConstraintCommon<3>
{
   public:
	using me_t = CConstraintMobileSlider;

	// Point indices: 0=mobile at the slider, 1 & 2: define the line on which 0
	// slides.

	CConstraintMobileSlider(
		const size_t _point_index, const size_t _ref_pt0, const size_t _ref_pt1)
		: CConstraintCommon({_point_index, _ref_pt0, _ref_pt1})
	{
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }
};

}  // namespace mbse
