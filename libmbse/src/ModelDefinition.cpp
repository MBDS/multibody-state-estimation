/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/ModelDefinition.h>
#include <mbse/AssembledRigidModel.h>
#include <mbse/constraints/ConstraintConstantDistance.h>
#include <mrpt/opengl.h>
#include <mrpt/core/round.h>
#include <mrpt/expr/CRuntimeCompiledExpression.h>

using namespace mbse;
using namespace std;

Body& ModelDefinition::addBody(const std::string& name)
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
	Body& new_body = *bodies_.rbegin();
	new_body.name = nam;

	return new_body;
}

/** Completely erases all defined points, joints, bodies, parameters, etc of
 * this object and leaves it blank. */
void ModelDefinition::clear() { *this = ModelDefinition(); }

MRPT_TODO("Initial position problem should refine these positions if needed.")

void ModelDefinition::setPointCoords(
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
void ModelDefinition::assembleRigidMBS(TSymbolicAssembledModel& armi) const
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
			const Body& b = bodies_[i];
			const auto nPts = b.points.size();

			switch (nPts)
			{
				case 0:
				case 1:
					THROW_EXCEPTION_FMT(
						"Body has an invalid number of points (=%u), valid are "
						">=2",
						static_cast<unsigned int>(nPts));
					break;
				case 2:
					ASSERT_(b.points[0] < points_.size());
					ASSERT_(b.points[1] < points_.size());
					const_cast<ModelDefinition*>(this)
						->addConstraint<ConstraintConstantDistance>(
							b.points[0], b.points[1], b.length());
					break;

				case 3:
				{
					const std::vector<std::pair<size_t, size_t>> pairs = {
						{0, 1}, {0, 2}, {1, 2}};

					for (const auto& p : pairs)
					{
						const size_t i = p.first, j = p.second;

						ASSERT_(b.points[i] < points_.size());
						ASSERT_(b.points[j] < points_.size());

						const double L = (points_.at(b.points[i]).coords -
										  points_.at(b.points[j]).coords)
											 .norm();

						const_cast<ModelDefinition*>(this)
							->addConstraint<ConstraintConstantDistance>(
								b.points[i], b.points[j], L);
					}
				}
				break;

				default:
					ASSERT_GE_(nPts, 4);
					THROW_EXCEPTION("TODO");
					break;
			};
		}
		// Mark these constraints as added:
		already_added_fixed_len_constraints_ = true;
	}
}

// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
std::shared_ptr<AssembledRigidModel> ModelDefinition::assembleRigidMBS(
	mrpt::optional_ref<const std::vector<RelativeDOF>> relativeCoordinates)
	const
{
	// 1) Build "symbolic" assembly:
	TSymbolicAssembledModel armi(*this);
	this->assembleRigidMBS(armi);

	// Append optional relative coordinates:
	if (relativeCoordinates.has_value()) armi.rDOFs = *relativeCoordinates;

	// 2) Actual assembly:
	return std::make_shared<AssembledRigidModel>(armi);
}

// -------------------------------------------------------------------
//                      FromYAML
// -------------------------------------------------------------------
ModelDefinition ModelDefinition::FromYAML(const mrpt::containers::yaml& c)
{
	using mrpt::expr::CRuntimeCompiledExpression;

	MRPT_START

	ASSERTMSG_(
		c.isMap(), mrpt::format(
					   "YAML node must be a map, but found: %s",
					   c.node().typeName().c_str()));

	ModelDefinition m;

	// ---------------------
	// Points
	// ---------------------
	ASSERT_(c["points"].isSequence());
	const auto yamlPts = c["points"];
	const auto nPts = yamlPts.asSequence().size();
	ASSERT_(nPts >= 1);
	m.setPointCount(nPts);

	CRuntimeCompiledExpression expX, expY;

	for (size_t idxPt = 0; idxPt < nPts; idxPt++)
	{
		const auto yamlPt = yamlPts.asSequence().at(idxPt).asMap();

		bool isFixed = false;
		if (auto it = yamlPt.find("fixed"); it != yamlPt.end())
			isFixed = (*it).second.as<bool>();

		const auto sX = yamlPt.at("x").as<std::string>();
		const auto sY = yamlPt.at("y").as<std::string>();

		expX.compile(sX);
		expY.compile(sY);

		m.setPointCoords(idxPt, {expX.eval(), expY.eval()}, isFixed);
	}

	// ---------------------
	// Planar bodies
	// ---------------------
	ASSERT_(c["planar_bodies"].isSequence());
	const auto yamlBodies = c["planar_bodies"];
	ASSERT_(yamlBodies.asSequence().size() >= 1);

	for (const auto& yamlBody : yamlBodies.asSequence())
	{
		CRuntimeCompiledExpression e;
		Body& b = m.addBody();

		std::map<std::string, double> expVars;
		expVars["index"] = m.bodies().size();

		MRPT_TODO("allow using x0,y0,x1,y1 as expressions.");

		const auto& yb = yamlBody.asMap();
		const auto pts = yb.at("points").asSequence();

		ASSERT_GE_(pts.size(), 2U);
		b.points.resize(pts.size());
		for (size_t ptIdx = 0; ptIdx < pts.size(); ptIdx++)
		{
			e.compile(
				pts.at(ptIdx).as<std::string>(), expVars,
				mrpt::format("points[%u]", static_cast<unsigned int>(ptIdx)));
			b.points.at(ptIdx) = mrpt::round(e.eval());
			ASSERT_LT_(b.points.at(ptIdx), m.getPointCount());
		}

		// Provide "auto" length:
		expVars["auto"] = (m.getPointInfo(b.points.at(0)).coords -
						   m.getPointInfo(b.points.at(1)).coords)
							  .norm();

		e.compile(yb.at("length").as<std::string>(), expVars, "length");
		b.length() = e.eval();

		expVars.erase("auto");
		expVars["length"] = b.length();

		e.compile(yb.at("mass").as<std::string>(), expVars, "mass");
		b.mass() = e.eval();
		expVars["mass"] = b.mass();

		e.compile(yb.at("I0").as<std::string>(), expVars, "I0");
		b.I0() = e.eval();
		expVars["I0"] = b.I0();

		ASSERT_EQUAL_(yb.at("cog").asSequence().size(), 2U);
		e.compile(
			yb.at("cog").asSequence().at(0).as<std::string>(), expVars,
			"cog.x");
		b.cog().x = e.eval();
		e.compile(
			yb.at("cog").asSequence().at(1).as<std::string>(), expVars,
			"cog.y");
		b.cog().y = e.eval();
	}

	// Constraints:
	// ---------------------
	MRPT_TODO("continue...");

	return m;

	MRPT_END
}
