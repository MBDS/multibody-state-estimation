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

#include "mbse-common.h"
#include <mrpt/opengl/CRenderizable.h>

namespace mbse
{
using namespace Eigen;
using namespace mrpt::math;

/** 2D generic body */
struct CBody
{
	CBody();  //!< Default ctor.

	std::string name;

	/** A 2D body is defined (in natural coords) with 2 points
	 * Indices of the body 2 points (from the list of all
	 * points in the problem); may include one fixed point
	 * (not a variable) */
	size_t points[2];

	/** In (kg) */
	inline double mass() const { return m_mass; }
	inline double& mass()
	{
		m_mass_matrices_cached = false;
		return m_mass;
	}

	/** Center of gravity (in local coordinates, origin=first point) */
	inline TPoint2D cog() const { return m_cog; }
	inline TPoint2D& cog()
	{
		m_mass_matrices_cached = false;
		return m_cog;
	}

	/** Fixed length (distance) between points 0-1 (constant since this is a
	 * rigid body) */
	inline double length() const { return m_length; }
	inline double& length()
	{
		m_mass_matrices_cached = false;
		return m_length;
	}

	/** Moment of inertia wrt point 0 */
	inline double I0() const { return m_I0; }
	inline double& I0()
	{
		m_mass_matrices_cached = false;
		return m_I0;
	}

	/** Computes the 3 different 2x2 blocks of the 4x4 mass matrix of a generic
	 * planar rigid element: \code [ M00   |  M01  ] M = [ ------+------ ] [
	 * M01^t |  M11  ] \endcode
	 */
	void evaluateMassMatrix(Matrix2d& M00, Matrix2d& M11, Matrix2d& M01) const;

	const Matrix2d& getM00() const;  //!< Computes (or gets the cached versions)
									 //!< of the mass submatrices
	const Matrix2d& getM11() const;  //!< Computes (or gets the cached versions)
									 //!< of the mass submatrices
	const Matrix2d& getM01() const;  //!< Computes (or gets the cached versions)
									 //!< of the mass submatrices

   private:
	/** Cached versions of mass submatrices, stored here after calling
	 * evaluateMassMatrix() */
	mutable Matrix2d m_M00, m_M11, m_M01;
	mutable bool m_mass_matrices_cached;

	/** Computes the 3 different 2x2 blocks of the 4x4 mass matrix of a generic
	 * planar rigid element */
	void internal_update_mass_submatrices() const;

	double m_mass;  //!< In (kg)
	TPoint2D m_cog;  //!< Center of gravity (in local coordinates, origin=first
					 //!< point)
	double m_length;  //!< Fixed length (distance) between points 0-1 (constant
					  //!< since this is a rigid body)
	double m_I0;  //!< Moment of inertia wrt point 0

   public:
	/**  Creates a 3D representation of the body */
	mrpt::opengl::CRenderizable::Ptr get3DRepresentation() const;

	/** Type of 3D object in which the body will be converted */
	enum render_style_t
	{
		reLine = 0,  //<! A simple line
		reCylinder  //<! A cylinder
	};

	struct TRenderParams
	{
		TRenderParams();  // Default ctor

		render_style_t render_style;  //!< Kind of object

		/* ==== Common options  ==== */
		bool show_grounds;  //!< Draws ground points as independent "ground
							//!< solids"
		double z_layer;  //!< Emulates links in "layers": an increment to be
						 //!< added to the Z coordinate of the object.

		/* ==== Render as lines ==== */
		uint8_t line_alpha;  //!< Transparency (0x00 - 0xff)
		float line_width;  //!< Line width (in pixels)

		/* ====  Render as cylinder ====  */
		double cyl_diameter;
	};

	TRenderParams render_params;

	// EIGEN_MAKE_ALIGNED_OPERATOR_NEW    // Required for aligned mem allocator
	// (only needed in classes containing fixed-size Eigen matrices)
};

}  // namespace mbse
