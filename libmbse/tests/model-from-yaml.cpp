/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <gtest/gtest.h>

#include <mbse/ModelDefinition.h>

TEST(ModelFromYaml, Pendulum)
{
	using namespace mrpt;

	const std::string sDef = R"(# Pendulum
points:
  - x: 0
    y: 0
    fixed: true
  - { x: 2.0*cos(deg2rad(35.0)), y: 2.0*sin(deg2rad(35.0)) }
planar_bodies:
  - points: [0, 1]
    length: auto
    mass: 1.0
    I0: (1/3)*mass*length^2
    cog: [0.5*length, 0.0]
)";

	const auto def = mrpt::containers::yaml::FromText(sDef);

	std::cout << "Model:\n" << def << std::endl;

	const auto model = mbse::ModelDefinition::FromYAML(def);

	EXPECT_EQ(model.getPointCount(), 2U);
	EXPECT_EQ(model.bodies().size(), 1U);

	EXPECT_EQ(model.getPointInfo(0).coords, mrpt::math::TPoint2D(0, 0));
	EXPECT_NEAR(model.getPointInfo(1).coords.x, 2 * cos(35.0_deg), 1e-5);
	EXPECT_NEAR(model.getPointInfo(1).coords.y, 2 * sin(35.0_deg), 1e-5);

	EXPECT_NEAR(model.bodies().at(0).length(), 2.0, 1e-5);
}
