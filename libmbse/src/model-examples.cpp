#include <mbse/model-examples.h>
#include <mrpt/random.h>

using namespace mbse;
using namespace mrpt::random;

void mbse::buildParameterizedMBS(
	const size_t nx, const size_t ny, CModelDefinition& model,
	const double NOISE_LEN)
{
	ASSERT_(nx >= 1 && ny >= 1);

#if MRPT_VERSION >= 0x199
	auto randomGenerator = mrpt::random::getRandomGenerator();
#else
	auto& randomGenerator = mrpt::random::randomGenerator;
#endif

	// Definition of constants
	const size_t npoints = (nx + 1) * (ny + 1);	 // number of points

	const double Lx = 2.0;	// Rod lengths
	const double Ly = 2.5;	// Rod lengths

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
			CBody& b = model.addBody();
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
	for (size_t row = 1; row <= ny;
		 row++)	 // from row 1 (0 belongs to fixed points) to ny
	{
		for (size_t j = 0; j <= nx; j++)  // from column 0 to nx
		{
			CBody& b = model.addBody();
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

	// cout << "# bodies: " << model.getBodies().size() << endl;
}

void mbse::buildLongStringMBS(const size_t N, CModelDefinition& model)
{
	ASSERT_(N >= 1);

	// Definition of constants
	const double L = 0.5;  // Rod lengths

	model.setPointCount(N + 1);

	// Definition of fixed points
	for (size_t i = 0; i <= N; i++)
		model.setPointCoords(i, TPoint2D(i * L, 0), i == 0 /*is fixed*/);

	// bars
	for (size_t j = 0; j < N; j++)
	{
		CBody& b = model.addBody();
		b.points[0] = j;
		b.points[1] = j + 1;
		b.length() = L;
		b.mass() = L * 0.1;
		b.I0() = b.mass() * mrpt::square(b.length()) / 3.0;
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
}

void mbse::buildFourBarsMBS(CModelDefinition& model)
{
	model.setPointCount(4);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 0));
	model.setPointCoords(2, TPoint2D(1, 2));
	model.setPointCoords(3, TPoint2D(4, 0), true /*is fixed*/);

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 1;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = 2;
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = -0.05;
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = std::sqrt(2.0 * 2.0 + 3.0 * 3.0);
		b.mass() = 4;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);

		b.render_params.z_layer = 0;
	}
}

void mbse::buildSliderCrankMBS(CModelDefinition& model)
{
	model.setPointCount(3);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(5, 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 1;
		b.points[1] = 2;

		b.length() = std::sqrt(17);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintFixedSlider(
		2 /*pt index*/, TPoint2D(-3, -2),
		TPoint2D(8, 2) /* The line on which to fix the point */
		));
}

void mbse::buildFollowerMBS(CModelDefinition& model)
{
	model.setPointCount(5);
	model.setPointCoords(0, TPoint2D(0, 0), true /*is fixed*/);
	model.setPointCoords(1, TPoint2D(1, 1));
	model.setPointCoords(2, TPoint2D(1, -4), true /*is fixed*/);
	model.setPointCoords(3, TPoint2D(1, 4));
	model.setPointCoords(4, TPoint2D(5, 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = std::sqrt(2);
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 2;
		b.points[1] = 3;

		b.length() = 8;
		b.mass() = 3;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}
	{
		CBody& b = model.addBody();
		b.points[0] = 3;
		b.points[1] = 4;

		b.length() = std::sqrt(32);
		b.mass() = 2;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintMobileSlider(
		1 /*pt index*/, 2, 3 /* two pts for defining the constraint */
		));

	model.addConstraint(CConstraintFixedSlider(
		4 /*pt index*/, TPoint2D(-5, 0),
		TPoint2D(10, 0) /* The line on which to fix the point */
		));
}

void mbse::buildTwoSliderBlocks(CModelDefinition& model)
{
	model.setPointCount(2);
	model.setPointCoords(0, TPoint2D(0, 15 * sin(mrpt::DEG2RAD(35))));
	model.setPointCoords(1, TPoint2D(15 * cos(mrpt::DEG2RAD(35)), 0));

	{
		CBody& b = model.addBody();
		b.points[0] = 0;
		b.points[1] = 1;

		b.length() = 15;
		b.mass() = 1;
		b.I0() = (1. / 3.) * b.mass() * mrpt::square(b.length());
		b.cog() = TPoint2D(b.length() * 0.5, 0);
	}

	model.addConstraint(CConstraintFixedSlider(
		0 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(0, 1) /* The line on which to fix the point */
		));
	model.addConstraint(CConstraintFixedSlider(
		1 /*pt index*/, TPoint2D(0, 0),
		TPoint2D(1, 0) /* The line on which to fix the point */
		));
}
