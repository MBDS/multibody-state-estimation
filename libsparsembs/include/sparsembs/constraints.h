#pragma once

#include "sparsembs-common.h"
#include <mrpt/opengl/CRenderizable.h>

namespace sparsembs
{
class CAssembledRigidModel;
using namespace Eigen;
using namespace mrpt::math;

/** The virtual base class of all constraint types. */
class CConstraintBase
{
   public:
	/** Alloc space for the needed rows, columns, etc. in the constraint vector
	 * and its Jacobians Each class should save a reference/pointer to the newly
	 * created elements, so they can be quickly updated in the \a update()
	 * method when invoked in the future. This method is called after the MBS is
	 * completely defined by the user, and before starting any kinematic/dynamic
	 * simulation.
	 */
	virtual void buildSparseStructures(CAssembledRigidModel& arm) const = 0;

	/** Update the previously allocated values given the current state of the
	 * MBS. This is called a very large number of times during simulations. */
	virtual void update(CAssembledRigidModel& arm) const = 0;

	/** Virtual destructor (required in any virtual base) */
	virtual ~CConstraintBase();

	virtual CConstraintBase* clone()
		const = 0;  //!< Clone operator for smart pointers

	/** Creates a 3D representation of the constraint, if applicable (e.g. the
	 * line of a fixed slider) \return false if the constraint has no 3D
	 * representation
	 */
	virtual bool get3DRepresentation(
		mrpt::opengl::CRenderizablePtr& inout_obj) const
	{
		return false;
	}
};

/** A smart pointer type for constraints */
typedef stlplus::smart_ptr_clone<CConstraintBase> CConstraintBasePtr;

/** Constraint: constant distance between two points */
class CConstraintConstantDistance : public CConstraintBase
{
   public:
	size_t point_index0, point_index1;
	double length;

	CConstraintConstantDistance(
		const size_t _point_index0, const size_t _point_index1,
		const double _length)
		: point_index0(_point_index0),
		  point_index1(_point_index1),
		  length(_length),
		  m_idx_constr(static_cast<size_t>(-1)),
		  dPhi_dx0(NULL),
		  dPhi_dy0(NULL),
		  dPhi_dx1(NULL),
		  dPhi_dy1(NULL),
		  dot_dPhi_dx0(NULL),
		  dot_dPhi_dy0(NULL),
		  dot_dPhi_dx1(NULL),
		  dot_dPhi_dy1(NULL)
	{
		m_points[0] = m_points[1] = NULL;
	}

	virtual void buildSparseStructures(CAssembledRigidModel& arm) const;
	virtual void update(CAssembledRigidModel& arm) const;

	virtual CConstraintBase* clone() const
	{
		return new CConstraintConstantDistance(
			point_index0, point_index1, length);
	}

   protected:
	mutable const TMBSPoint* m_points[2];
	mutable TPoint2DOF
		m_pointDOFs[2];  // The indices of each point in the state vector "q"
	mutable size_t m_idx_constr;
	mutable double *dPhi_dx0, *dPhi_dy0, *dPhi_dx1,
		*dPhi_dy1;  // Pointers to entries in the sparse Jacobian dPhi_dq
	mutable double *dot_dPhi_dx0, *dot_dPhi_dy0, *dot_dPhi_dx1,
		*dot_dPhi_dy1;  // Pointers to entries in the sparse Jacobian
						// \dot{dPhi_dq}
};

/** Constraint: forces a point to lie exactly on a fixed line (e.g. sliders) */
class CConstraintFixedSlider : public CConstraintBase
{
   public:
	size_t point_index;
	TPoint2D line_pt[2];  //!< The point is forced to lie on the line defined by
						  //!< these two fixed points (x0,y0)-(x1,y1)

	CConstraintFixedSlider(
		const size_t _point_index, const TPoint2D& pt0, const TPoint2D& pt1)
		: point_index(_point_index),
		  m_point(NULL),
		  m_idx_constr(static_cast<size_t>(-1)),
		  dPhi_dx0(NULL),
		  dPhi_dy0(NULL),
		  dot_dPhi_dx0(NULL),
		  dot_dPhi_dy0(NULL)
	{
		line_pt[0] = pt0;
		line_pt[1] = pt1;
	}

	virtual void buildSparseStructures(CAssembledRigidModel& arm) const;
	virtual void update(CAssembledRigidModel& arm) const;

	virtual CConstraintBase* clone() const
	{
		return new CConstraintFixedSlider(point_index, line_pt[0], line_pt[1]);
	}

	/** Creates a 3D representation of the constraint, if applicable (e.g. the
	 * line of a fixed slider) \return false if the constraint has no 3D
	 * representation
	 */
	virtual bool get3DRepresentation(
		mrpt::opengl::CRenderizablePtr& inout_obj) const;

   protected:
	mutable const TMBSPoint* m_point;
	mutable TPoint2DOF
		m_pointDOF;  // The indices of each point in the state vector "q"
	mutable size_t m_idx_constr;
	mutable TPoint2D m_Delta;  //!< pt[1]-pt[0], precomputed only once
	mutable double *dPhi_dx0,
		*dPhi_dy0;  // Pointers to entries in the sparse Jacobian dPhi_dq
	mutable double *dot_dPhi_dx0,
		*dot_dPhi_dy0;  // Pointers to entries in the sparse Jacobian
						// \dot{dPhi_dq}
};

/** Constraint: forces a point to lie exactly on the line defined by two other
 * reference points */
class CConstraintMobileSlider : public CConstraintBase
{
   public:
	size_t point_index;
	size_t ref_pts[2];  //!< The point is forced to lie on the line defined by
						//!< these two points

	CConstraintMobileSlider(
		const size_t _point_index, const size_t _ref_pt0, const size_t _ref_pt1)
		: point_index(_point_index), m_idx_constr(static_cast<size_t>(-1))
	{
		ref_pts[0] = _ref_pt0;
		ref_pts[1] = _ref_pt1;
		for (int i = 0; i < 3; i++)
		{
			m_points[i] = NULL;
			dPhi_dx[i] = dPhi_dy[i] = NULL;
			dot_dPhi_dx[i] = dot_dPhi_dy[i] = NULL;
		}
	}

	virtual void buildSparseStructures(CAssembledRigidModel& arm) const;
	virtual void update(CAssembledRigidModel& arm) const;

	virtual CConstraintBase* clone() const
	{
		return new CConstraintMobileSlider(point_index, ref_pts[0], ref_pts[1]);
	}

   protected:
	mutable const TMBSPoint* m_points[3];
	mutable TPoint2DOF
		m_pointDOF[3];  // The indices of each point in the state vector "q"
	mutable size_t m_idx_constr;
	mutable double *dPhi_dx[3],
		*dPhi_dy[3];  // Pointers to entries in the sparse Jacobian dPhi_dq
	mutable double *dot_dPhi_dx[3],
		*dot_dPhi_dy[3];  // Pointers to entries in the sparse Jacobian
						  // \dot{dPhi_dq}
};

}  // namespace sparsembs
