/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2021 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/model-examples.h>
#include <mbse/constraints/ConstraintFixedSlider.h>
#include <mbse/constraints/ConstraintMobileSlider.h>
#include <mrpt/random.h>

using namespace mbse;
using namespace mrpt::random;
using namespace mrpt::math;

ModelDefinition mbse::buildParameterizedMBS(
	const size_t nx, const size_t ny, const double NOISE_LEN)
{
	ASSERT_(nx >= 1 && ny >= 1);

	auto randomGenerator = mrpt::random::getRandomGenerator();

	// Definition of constants
	const size_t npoints = (nx + 1) * (ny + 1);	 // number of points

	const double Lx = 2.0;	// Rod lengths
	const double Ly = 2.5;	// Rod lengths

	ModelDefinition model;
	model.setPointCount(npoints);

	// Definition of fixed points
	for (size_t i = 0; i <= nx; i++)
	{
		const double x =
			i * Lx;	 // + (irregular_mesh ?
					 // randomGenerator.drawUniform(-NOISE_LEN,NOISE_LEN) : 0 );
		model.setPointCoords(i, TPoint2D(x, 0), true /*is fixed*/);
	}

	// points definition
	for (size_t row = 1; row <= ny;
		 row++)	 // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j <= nx; j++)  // from column 0 to nx
		{
			const double x =
				j * Lx + randomGenerator.drawUniform(-NOISE_LEN, NOISE_LEN);
			const double y =
				row * Ly + randomGenerator.drawUniform(-NOISE_LEN, NOISE_LEN);
			model.setPointCoords(row * (nx + 1) + j, TPoint2D(x, y));
		}
	}

	// horizontal bars
	for (size_t row = 1; row <= ny;
		 row++)	 // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j < nx; j++)	 // from column 0 to nx-1
		{
			Body& b = model.addBody();
			b.points[0] = row * (nx + 1) + j;  // left side point
			b.points[1] = row * (nx + 1) + j + 1;  // right side point
			b.length() = (model.getPointInfo(b.points[0]).coords -
						  model.getPointInfo(b.points[1]).coords)
							 .norm();
			b.mass() = 1 * b.length();
			b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
			b.cog() = TPoint2D(b.length() * 0.5, 0);
		}
	}

	// vertical bars
	// from row 1 (0 belongs to fixed points) to ny
	for (size_t row = 1; row <= ny; row++)
	{
		for (size_t j = 0; j <= nx; j++)  // from column 0 to nx
		{
			Body& b = model.addBody();
			b.points[0] = row * (nx + 1) + j;  // upper point
			b.points[1] = row * (nx + 1) + j - nx - 1;	// lower point
			b.length() = (model.getPointInfo(b.points[0]).coords -
						  model.getPointInfo(b.points[1]).coords)
							 .norm();
			b.mass() = 1 * b.length();
			b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
			b.cog() = TPoint2D(b.length() * 0.5, 0);

			b.render_params.z_layer = 0.1;
		}
	}

	return model;
}

ModelDefinition mbse::buildLongStringMBS(
	const size_t N, double segmentLength, double segmentMassPerMeter)
{
	ASSERT_(N >= 1);

	// Definition of constants
	const double L = segmentLength;	 // Rod lengths

	ModelDefinition model;
	model.setPointCount(N + 1);

	// Definition of fixed points
	for (size_t i = 0; i <= N; i++)
		model.setPointCoords(i, TPoint2D(i * L, 0), i == 0 /*is fixed*/);

	// bars
	for (size_t j = 0; j < N; j++)
	{
		Body& b = model.addBody();
		b.points[0] = j;
		b.points[1] = j + 1;
		b.length() = L;
		b.mass() = L * segmentMassPerMeter;
		b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	return model;
}

ModelDefinition mbse::buildFourBarsMBS()
{
	ModelDefinition model;
	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 0));
	model.setPointCoords(2, TPoint2D(1, 2));
	model.setPointCoords(3, TPoint2D(4, 0), true /*is fixed*/);

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 1;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	{
		Body& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = -0.05;
	}
	{
		Body& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = std::sqrt(2.0 * 2.0 + 3.0 * 3.0);
		b.mass() = 4;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	return model;
}

ModelDefinition mbse::buildSliderCrankMBS()
{
	ModelDefinition model;
	model.setPointCount(3);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(5, 0));

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		Body& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = std::sqrt(17);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint<ConstraintFixedSlider>(
		2 /*pt index*/, TPoint2D(-3, -2),
		TPoint2D(8, 2) /* The line on which to fix the point */
	);
	return model;
}

ModelDefinition mbse::buildFollowerMBS()
{
	ModelDefinition model;
	model.setPointCount(5);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(1, -4), true /*is fixed*/);
	model.setPointCoords(3, TPoint2D(1, 4));
	model.setPointCoords(4, TPoint2D(5, 0));

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		Body& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 8;
		b.mass() = 3;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		Body& b = model.addBody();
		b.points[0] = 3;
		b.points[1] = 4;

		b.length() = std::sqrt(32);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint<ConstraintMobileSlider>(
		1 /*pt index*/, 2, 3 /* two pts for defining the constraint */
	);

	model.addConstraint<ConstraintFixedSlider>(
		4 /*pt index*/, TPoint2D(-5, 0),
		TPoint2D(10, 0) /* The line on which to fix the point */
	);
	return model;
}

ModelDefinition mbse::buildTwoSliderBlocks()
{
	ModelDefinition model;
	model.setPointCount(2);
	model.setPointCoords(0, TPoint2D(0, 15 * sin(mrpt::DEG2RAD(35))));
	model.setPointCoords(1, TPoint2D(15 * cos(mrpt::DEG2RAD(35)), 0));

	{
		Body& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 15;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint<ConstraintFixedSlider>(
		0 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(0, 1) /* The line on which to fix the point */
	);
	model.addConstraint<ConstraintFixedSlider>(
		1 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(1, 0) /* The line on which to fix the point */
	);
	return model;
}
