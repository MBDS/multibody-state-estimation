/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/ModelDefinition.h>
#include <mbse/AssembledRigidModel.h>
#include <mbse/constraints/ConstraintConstantDistance.h>
#include <mbse/constraints/ConstraintRelativePosition.h>
#include <mrpt/opengl.h>
#include <mrpt/core/round.h>
#include <mrpt/math/TPose2D.h>
#include <mrpt/expr/CRuntimeCompiledExpression.h>
#include <cmath>  // atan2()

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
				{
					THROW_EXCEPTION_FMT(
						"Body has an invalid number of points (=%u), valid are "
						">=2",
						static_cast<unsigned int>(nPts));
				}
				break;
				case 2:
				{
					ASSERT_(b.points[0] < points_.size());
					ASSERT_(b.points[1] < points_.size());
					const_cast<ModelDefinition*>(this)
						->addConstraint<ConstraintConstantDistance>(
							b.points[0], b.points[1], b.length());
				}
				break;

				case 3:
				{
					MRPT_TODO("Handle 3 aligned points");

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

#if 0
				default:
				{
					ASSERT_GE_(nPts, 4);

					// 3 first points as in case above:
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

					for (size_t i = 3; i < nPts; i++)
					{
						const_cast<ModelDefinition*>(this)
							->addConstraint<ConstraintRelativePosition>(
								b.points[0], b.points[1], b.points[2],
								b.points[i]);
					}

					break;
				}
#else
				default:
				{
					MRPT_TODO("Handle 3 aligned points");

					std::vector<std::pair<size_t, size_t>> pairs = {{0, 1}};

					for (size_t j = 2; j < nPts; j++)
					{
						pairs.emplace_back(0, j);
						pairs.emplace_back(1, j);
					}

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

#endif
			};
		}
		// Mark these constraints as added:
		already_added_fixed_len_constraints_ = true;
	}

	// Append optional relative coordinates:
	armi.rDOFs = rDOFs_;
}

// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
std::shared_ptr<AssembledRigidModel> ModelDefinition::assembleRigidMBS() const
{
	// 1) Build "symbolic" assembly:
	TSymbolicAssembledModel armi(*this);
	this->assembleRigidMBS(armi);

	// 2) Actual assembly:
	auto arm = std::make_shared<AssembledRigidModel>(armi);

	arm->realize_operating_point();

	return arm;
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
	// Parameters
	// ---------------------
	std::map<std::string, double> expVars;	// expression variables
	const std::set<std::string> reservedNames = {
		"index", "auto", "length", "mass", "I0"};

	if (c.has("parameters"))
	{
		ASSERT_(c["parameters"].isMap());
		for (const auto& kv : c["parameters"].asMap())
		{
			const auto name = kv.first.as<std::string>();
			const auto value = kv.second.as<std::string>();

			ASSERTMSG_(
				reservedNames.count(name) == 0,
				mrpt::format(
					"Cannot use reserved keyword '%s' as parameter name.",
					name.c_str()));

			CRuntimeCompiledExpression exp;
			exp.compile(value, expVars, name);

			expVars[name] = exp.eval();
		}
	}

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

		expX.compile(sX, expVars);
		expY.compile(sY, expVars);

		const auto pt = mrpt::math::TPoint2D(expX.eval(), expY.eval());

		m.setPointCoords(idxPt, pt, isFixed);

		expVars[std::string("x") + std::to_string(idxPt)] = pt.x;
		expVars[std::string("y") + std::to_string(idxPt)] = pt.y;
	}

	// ---------------------
	// Planar bodies
	// ---------------------
	ASSERT_(c["planar_bodies"].isSequence());
	const auto yamlBodies = c["planar_bodies"];
	ASSERT_(yamlBodies.asSequence().size() >= 1);

	for (const auto& yamlBody : yamlBodies.asSequence())
	{
		Body& b = m.addBody();

		expVars["index"] = m.bodies().size();

		const auto& yb = yamlBody.asMap();
		const auto pts = yb.at("points").asSequence();

		ASSERT_GE_(pts.size(), 2U);
		b.points.resize(pts.size());
		for (size_t ptIdx = 0; ptIdx < pts.size(); ptIdx++)
		{
			CRuntimeCompiledExpression e;
			e.compile(
				pts.at(ptIdx).as<std::string>(), expVars,
				mrpt::format("points[%u]", static_cast<unsigned int>(ptIdx)));
			b.points.at(ptIdx) = mrpt::round(e.eval());
			ASSERT_LT_(b.points.at(ptIdx), m.getPointCount());
		}

		// Provide "auto" length:
		// "length" is always "auto" for N>=3 bodies:
		if (pts.size() < 3)
		{
			expVars["auto"] = (m.getPointInfo(b.points.at(0)).coords -
							   m.getPointInfo(b.points.at(1)).coords)
								  .norm();

			CRuntimeCompiledExpression e;
			e.compile(yb.at("length").as<std::string>(), expVars, "length");
			b.length() = e.eval();

			expVars.erase("auto");
			expVars["length"] = b.length();
		}
		else
		{
			// 3-ary or higher-order: provide helper length variables:
			for (unsigned int i = 0; i < pts.size(); i++)
			{
				for (unsigned int j = i + 1; j < pts.size(); j++)
				{
					const double L = (m.getPointInfo(b.points.at(i)).coords -
									  m.getPointInfo(b.points.at(j)).coords)
										 .norm();
					expVars[mrpt::format("length%u%u", i + 1, j + 1)] = L;
					expVars[mrpt::format("length%u%u", j + 1, i + 1)] = L;

					// We need L01 for some fixed formulas later on:
					if (i == 0 && j == 1) b.length() = L;
				}
			}
		}
		ASSERT_(b.length() > 0);

		// Build local point coordinates:
		{
			auto& localPts = b.fixedPointsLocal();

			const auto& pt0 = m.getPointInfo(b.points.at(0)).coords;
			const auto& pt1 = m.getPointInfo(b.points.at(1)).coords;
			const auto v01 = pt1 - pt0;
			ASSERT_(v01.norm() > 0);

			const auto refPose =
				mrpt::math::TPose2D(pt0.x, pt0.y, std::atan2(v01.y, v01.x));

			localPts.clear();
			for (size_t i = 0; i < pts.size(); i++)
			{
				const auto localPt = refPose.inverseComposePoint(
					m.getPointInfo(b.points.at(i)).coords);
				localPts.push_back(localPt);
			}
		}

		CRuntimeCompiledExpression e;
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

		expVars.erase("mass");
		expVars.erase("length");
		expVars.erase("I0");
	}

	// relative coordinates
	// ------------------------------------
	if (c.has("relative_coordinates"))
	{
		ASSERT_(c["relative_coordinates"].isSequence());
		for (const auto& rc : c["relative_coordinates"].asSequence())
		{
			ASSERT_(rc.isMap());
			const auto& rcm = rc.asMap();
			ASSERT_(rcm.count("type"));
			ASSERT_(rcm.count("points"));
			ASSERT_(rcm.at("points").isSequence());

			// parse:
			const auto type = rcm.at("type").as<std::string>();
			std::vector<size_t> ptIndices;
			for (const auto& e : rcm.at("points").asSequence())
				ptIndices.push_back(e.as<size_t>());

			// create rDoF:
			if (type == "RelativeAngleAbsoluteDOF")
			{
				ASSERT_EQUAL_(ptIndices.size(), 2);
				m.rDOFs_.emplace_back(
					RelativeAngleAbsoluteDOF(ptIndices[0], ptIndices[1]));
			}
			else
			{
				THROW_EXCEPTION_FMT(
					"Unknown relative coordinate type: '%s'", type.c_str());
			}
		}
	}

	// Constraints:
	// ---------------------
	MRPT_TODO("continue...");

	return m;

	MRPT_END
}
