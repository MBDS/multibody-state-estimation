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
/** Constraint: forces a point to lie exactly on a fixed line (e.g. sliders) */
class CConstraintFixedSlider : public CConstraintBase,
							   public CConstraintCommon<1>
{
   public:
	using me_t = CConstraintFixedSlider;

	mrpt::math::TPoint2D
		line_pt[2];  //!< The point is forced to lie on the line defined by
					 //!< these two fixed points (x0,y0)-(x1,y1)

	CConstraintFixedSlider(
		const size_t _point_index, const mrpt::math::TPoint2D& pt0,
		const mrpt::math::TPoint2D& pt1)
		: CConstraintCommon({_point_index})
	{
		line_pt[0] = pt0;
		line_pt[1] = pt1;
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	Ptr clone() const override { return std::make_shared<me_t>(*this); }

	bool get3DRepresentation(
		mrpt::opengl::CRenderizable::Ptr& inout_obj) const override;

   protected:
	/** pt[1]-pt[0], precomputed only once */
	mutable mrpt::math::TPoint2D Delta_;
};

}  // namespace mbse
