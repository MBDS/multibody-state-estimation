#include <sparsembs/CBody.h>
#include <mrpt/opengl.h>

using namespace sparsembs;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt::utils;

MRPT_TODO("Allow bodies with more than 2 points")

// Constructor:
CBody::CBody()
	: name(),
	  m_mass_matrices_cached(false),
	  m_mass(0),
	  m_cog(0, 0),
	  m_length(0),
	  m_I0(0)
{
	points[0] = points[1] = static_cast<size_t>(-1);
}

/* Computes the 3 different 2x2 blocks of the 4x4 mass matrix of a generic
 * planar rigid element: [ M00   |  M01  ] M = [ ------+------ ] [ M01^t |  M11
 * ] In fact it's enough to fill the upper-triangular parts of M00 and M11. M01
 * is not symmetric so it must be filled-in entirely.
 */
void CBody::evaluateMassMatrix(
	Matrix2d& M00, Matrix2d& M11, Matrix2d& M01) const
{
	if (!m_mass_matrices_cached) internal_update_mass_submatrices();

	M00 = m_M00;
	M01 = m_M01;
	M11 = m_M11;
}

/** Computes the 3 different 2x2 blocks of the 4x4 mass matrix of a generic
 * planar rigid element */
void CBody::internal_update_mass_submatrices() const
{
	// From ch3, eq(17) of "ANALISIS Y SINTESIS DE MECANISMOS POR ORDENADOR",
	// Javier Cuadrado
	const double MXg_L = m_mass * m_cog.x / m_length;
	const double MYg_L = m_mass * m_cog.y / m_length;
	const double I0_L2 = m_I0 / (m_length * m_length);

	// M00:
	m_M00(0, 0) = m_mass - 2 * MXg_L + I0_L2;
	m_M00(0, 1) = 0;
	m_M00(1, 0) = 0;
	m_M00(1, 1) = m_M00(0, 0);

	// M11:
	m_M11(0, 0) = I0_L2;
	m_M11(0, 1) = 0;
	m_M11(1, 0) = 0;
	m_M11(1, 1) = m_M11(0, 0);

	// M01:
	m_M01(0, 0) = MXg_L - I0_L2;
	m_M01(0, 1) = -MYg_L;
	m_M01(1, 0) = MYg_L;
	m_M01(1, 1) = m_M01(0, 0);

	m_mass_matrices_cached = true;
}

const Matrix2d& CBody::getM00() const
{
	if (!m_mass_matrices_cached) internal_update_mass_submatrices();
	return m_M00;
}
const Matrix2d& CBody::getM11() const
{
	if (!m_mass_matrices_cached) internal_update_mass_submatrices();
	return m_M11;
}
const Matrix2d& CBody::getM01() const
{
	if (!m_mass_matrices_cached) internal_update_mass_submatrices();
	return m_M01;
}

/**  Creates a 3D representation of the body */
mrpt::opengl::CRenderizablePtr CBody::get3DRepresentation() const
{
	mrpt::opengl::CSetOfObjectsPtr objs = mrpt::opengl::CSetOfObjects::Create();

	switch (render_params.render_style)
	{
		case reCylinder:
		{
			mrpt::opengl::CCylinderPtr obj = mrpt::opengl::CCylinder::Create();

			obj->setRadius(render_params.cyl_diameter);
			obj->setColor_u8(mrpt::utils::TColor(0xFF, 0x00, 0x00));

			obj->setHeight(this->m_length);

			// Cylinder aligned with +X
			obj->setPose(mrpt::poses::CPose3D(
				0, 0, render_params.z_layer, DEG2RAD(0), DEG2RAD(90),
				DEG2RAD(0)));

			objs->insert(obj);
		}
		break;

		case reLine:
		{
			mrpt::opengl::CSimpleLinePtr obj =
				mrpt::opengl::CSimpleLine::Create();

			obj->setLineWidth(render_params.line_width);
			obj->enableAntiAliasing(true);
			obj->setColor_u8(mrpt::utils::TColor(
				0xFF, 0xFF, 0xFF, render_params.line_alpha));

			obj->setLineCoords(
				0, 0, render_params.z_layer, this->m_length, 0,
				render_params.z_layer);

			objs->insert(obj);
		}
		break;

		default:
			THROW_EXCEPTION("Unhandled render style!?")
	}

	return objs;
}

CBody::TRenderParams::TRenderParams()
	: render_style(reCylinder),
	  show_grounds(true),
	  /* ==== Common options  ==== */
	  z_layer(0),
	  /* ==== Render as lines ==== */
	  line_alpha(0x8f),
	  line_width(0.8f),
	  /* ====  Render as cylinder ====  */
	  cyl_diameter(0.05)
{
}
