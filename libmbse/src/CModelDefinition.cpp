/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/CModelDefinition.h>
#include <mbse/CAssembledRigidModel.h>
#include <mrpt/opengl.h>

using namespace mbse;
using namespace std;

CModelDefinition::CModelDefinition()
	: already_added_fixed_len_constraints_(false)
{
}

CBody& CModelDefinition::addBody(const std::string& name)
{
	ASSERTMSG_(
		!already_added_fixed_len_constraints_,
		"Can't modify model after assembling!");

	// Build name if none provided:
	const std::string nam =
		name.empty()
			? mrpt::format("body%u", static_cast<unsigned int>(bodies_.size()))
			: name;

	// Create:
	bodies_.resize(bodies_.size() + 1);

	// Set name & return:
	CBody& new_body = *bodies_.rbegin();
	new_body.name = nam;

	return new_body;
}

/** Completely erases all defined points, joints, bodies, parameters, etc of
 * this object and leaves it blank. */
void CModelDefinition::clear() { *this = CModelDefinition(); }

MRPT_TODO("Initial position problem should refine these positions if needed.")

void CModelDefinition::setPointCoords(
	const size_t i, const mrpt::math::TPoint2D& coords, const bool is_fixed)
{
	ASSERT_(i < points_.size());

	Point2& pt = points_[i];
	pt.coords = coords;
	pt.fixed = is_fixed;
}

// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
void CModelDefinition::assembleRigidMBS(TSymbolicAssembledModel& armi) const
{
	armi.clear();

	// 1) Count number of natural coordinates which are unknowns ==> # of scalar
	// unknowns in vector q
	// ---------------------------------------------------------------------------------------------------
	const size_t nPts = points_.size();
	for (size_t i = 0; i < nPts; i++)
	{
		const Point2& pt = points_[i];
		// If it's not fixed, add to the list of DOFs
		if (!pt.fixed)
		{
			armi.DOFs.emplace_back(i, PointDOF::X);
			armi.DOFs.emplace_back(i, PointDOF::Y);
		}
	}

	// 2) For each rigid body, automatically add one constant-distance
	// constraint:
	// ---------------------------------------------------------------------------------------------------
	if (!already_added_fixed_len_constraints_)
	{
		for (size_t i = 0; i < bodies_.size(); i++)
		{
			const CBody& b = bodies_[i];
			ASSERT_(b.points[0] < points_.size());
			ASSERT_(b.points[1] < points_.size());

			const_cast<CModelDefinition*>(this)->addConstraint(
				CConstraintConstantDistance(
					b.points[0], b.points[1], b.length()));
		}
		// Mark these constraints as added:
		already_added_fixed_len_constraints_ = true;
	}
}

// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
std::shared_ptr<CAssembledRigidModel> CModelDefinition::assembleRigidMBS(
	mrpt::optional_ref<const std::vector<RelativeDOF>> relativeCoordinates)
{
	// 1) Build "symbolic" assembly:
	TSymbolicAssembledModel armi(*this);
	this->assembleRigidMBS(armi);

	// Append optional relative coordinates:
	if (relativeCoordinates.has_value()) armi.rDOFs = *relativeCoordinates;

	// 2) Actual assembly:
	return std::make_shared<CAssembledRigidModel>(armi);
}
