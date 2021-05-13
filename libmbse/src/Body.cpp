/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/Body.h>
#include <mrpt/opengl.h>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;

MRPT_TODO("Allow bodies with more than 2 points")

void Body::evaluateMassMatrix(Matrix2d& M00, Matrix2d& M11, Matrix2d& M01) const
{
	if (!mass_matrices_cached_) internal_update_mass_submatrices();

	M00 = M00_;
	M01 = M01_;
	M11 = M11_;
}

/** Computes the 3 different 2x2 blocks of the 4x4 mass matrix of a generic
 * planar rigid element */
void Body::internal_update_mass_submatrices() const
{
	// From ch3, eq(17) of "ANALISIS Y SINTESIS DE MECANISMOS POR ORDENADOR",
	// Javier Cuadrado
	const double MXg_L = mass_ * cog_.x / length_;
	const double MYg_L = mass_ * cog_.y / length_;
	const double I0_L2 = I0_ / (length_ * length_);

	// M00:
	M00_(0, 0) = mass_ - 2 * MXg_L + I0_L2;
	M00_(0, 1) = 0;
	M00_(1, 0) = 0;
	M00_(1, 1) = M00_(0, 0);

	// M11:
	M11_(0, 0) = I0_L2;
	M11_(0, 1) = 0;
	M11_(1, 0) = 0;
	M11_(1, 1) = M11_(0, 0);

	// M01:
	M01_(0, 0) = MXg_L - I0_L2;
	M01_(0, 1) = -MYg_L;
	M01_(1, 0) = MYg_L;
	M01_(1, 1) = M01_(0, 0);

	mass_matrices_cached_ = true;
}

const Matrix2d& Body::getM00() const
{
	if (!mass_matrices_cached_) internal_update_mass_submatrices();
	return M00_;
}
const Matrix2d& Body::getM11() const
{
	if (!mass_matrices_cached_) internal_update_mass_submatrices();
	return M11_;
}
const Matrix2d& Body::getM01() const
{
	if (!mass_matrices_cached_) internal_update_mass_submatrices();
	return M01_;
}

/**  Creates a 3D representation of the body */
mrpt::opengl::CRenderizable::Ptr Body::get3DRepresentation() const
{
	using namespace mrpt;

	mrpt::opengl::CSetOfObjects::Ptr objs =
		mrpt::opengl::CSetOfObjects::Create();

	switch (render_params.render_style)
	{
		case reCylinder:
		{
			auto obj = mrpt::opengl::CCylinder::Create();

			obj->setRadius(render_params.cyl_diameter);
			obj->setColor_u8(mrpt::img::TColor(0xFF, 0x00, 0x00));

			obj->setHeight(this->length_);

			// Cylinder aligned with +X
			obj->setPose(mrpt::poses::CPose3D(
				0, 0, render_params.z_layer, DEG2RAD(0), DEG2RAD(90),
				DEG2RAD(0)));

			objs->insert(obj);
		}
		break;

		case reLine:  // synonum with: reSimplex:
		{
			ASSERT_EQUAL_(fixedPointsLocal_.size(), points.size());

			if (points.size() == 2)
			{
				auto obj = mrpt::opengl::CSimpleLine::Create();

				obj->setLineWidth(render_params.line_width);
				obj->enableAntiAliasing(true);
				obj->setColor_u8(mrpt::img::TColor(
					0xFF, 0xFF, 0xFF, render_params.line_alpha));

				obj->setLineCoords(
					0, 0, render_params.z_layer, this->length_, 0,
					render_params.z_layer);

				objs->insert(obj);
			}
			else
			{
				auto obj = mrpt::opengl::CSetOfLines::Create();

				obj->setLineWidth(render_params.line_width);
				obj->enableAntiAliasing(true);
				obj->setColor_u8(mrpt::img::TColor(
					0xFF, 0xFF, 0xFF, render_params.line_alpha));

				obj->appendLine(
					fixedPointsLocal_.at(0).x, fixedPointsLocal_.at(0).y,
					render_params.z_layer,	//
					fixedPointsLocal_.at(1).x, fixedPointsLocal_.at(1).y,
					render_params.z_layer);
				for (size_t i = 2; i <= points.size(); i++)
				{
					const size_t j = i % points.size();
					obj->appendLineStrip(
						fixedPointsLocal_.at(j).x, fixedPointsLocal_.at(j).y,
						render_params.z_layer);
				}

				objs->insert(obj);
			}
		}
		break;

		default:
			THROW_EXCEPTION("Unhandled render style!?");
	}

	return objs;
}
