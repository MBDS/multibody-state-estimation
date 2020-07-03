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

#include <mbse/mbse-common.h>
#include <mrpt/opengl/CRenderizable.h>

namespace mbse
{
class CAssembledRigidModel;

/** The virtual base class of all constraint types. */
class CConstraintBase
{
   public:
	/** A smart pointer type for constraints */
	using Ptr = std::shared_ptr<CConstraintBase>;

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

	/** Clone operator for smart pointers */
	virtual CConstraintBase* clone() const = 0;

	/** Creates a 3D representation of the constraint, if applicable (e.g.
   the
	 * line of a fixed slider) \return false if the constraint has no 3D
	 * representation
	 */
	virtual bool get3DRepresentation(
		[[maybe_unused]]  //
		mrpt::opengl::CRenderizable::Ptr& inout_obj) const
	{
		return false;
	}
};

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
		  length(_length)
	{
	}

	virtual void buildSparseStructures(CAssembledRigidModel& arm) const;
	virtual void update(CAssembledRigidModel& arm) const;

	virtual CConstraintBase* clone() const
	{
		return new CConstraintConstantDistance(*this);
	}

   protected:
	mutable const Point2* m_points[2] = {nullptr, nullptr};
	/** The indices of each point in the state vector "q" */
	mutable Point2ToDOF m_pointDOFs[2];
	mutable size_t m_idx_constr = static_cast<size_t>(-1);

	/** Pointers to entries in the sparse Jacobian dPhi_dq */
	mutable double *dPhi_dx0 = nullptr, *dPhi_dy0 = nullptr,
				   *dPhi_dx1 = nullptr, *dPhi_dy1 = nullptr;

	/** Pointers to entries in the sparse Jacobian \dot{dPhi_dq} */
	mutable double *dot_dPhi_dx0 = nullptr, *dot_dPhi_dy0 = nullptr,
				   *dot_dPhi_dx1 = nullptr, *dot_dPhi_dy1 = nullptr;

	/** Pointers to entries in the sparse Jacobian d(Phiq*dq)_dq */
	mutable double *dPhiqdq_dx0 = nullptr, *dPhiqdq_dy0 = nullptr,
				   *dPhiqdq_dx1 = nullptr, *dPhiqdq_dy1 = nullptr;
};

/** Constraint: forces a point to lie exactly on a fixed line (e.g. sliders) */
class CConstraintFixedSlider : public CConstraintBase
{
   public:
	size_t point_index;
	mrpt::math::TPoint2D
		line_pt[2];  //!< The point is forced to lie on the line defined by
					 //!< these two fixed points (x0,y0)-(x1,y1)

	CConstraintFixedSlider(
		const size_t _point_index, const mrpt::math::TPoint2D& pt0,
		const mrpt::math::TPoint2D& pt1)
		: point_index(_point_index)
	{
		line_pt[0] = pt0;
		line_pt[1] = pt1;
	}

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	CConstraintBase* clone() const override
	{
		return new CConstraintFixedSlider(*this);
	}

	bool get3DRepresentation(
		mrpt::opengl::CRenderizable::Ptr& inout_obj) const override;

   protected:
	mutable const Point2* m_point = nullptr;
	// The indices of each point in the state vector "q"
	mutable Point2ToDOF m_pointDOF;
	mutable size_t m_idx_constr = static_cast<size_t>(-1);
	mutable mrpt::math::TPoint2D
		m_Delta;  //!< pt[1]-pt[0], precomputed only once

	/** Pointers to entries in the sparse Jacobian dPhi_dq */
	mutable double *dPhi_dx0 = nullptr, *dPhi_dy0 = nullptr;

	/** Pointers to entries in the sparse Jacobian \dot{dPhi_dq} */
	mutable double *dot_dPhi_dx0 = nullptr, *dot_dPhi_dy0 = nullptr;

	/** Pointers to entries in the sparse Jacobian d(Phiq*dq)_dq */
	mutable double *dPhiqdq_dx0 = nullptr, *dPhiqdq_dy0 = nullptr;
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

	void buildSparseStructures(CAssembledRigidModel& arm) const override;
	void update(CAssembledRigidModel& arm) const override;

	CConstraintBase* clone() const override
	{
		return new CConstraintMobileSlider(*this);
	}

   protected:
	mutable const Point2* m_points[3];
	// The indices of each point in the state vector "q"
	mutable Point2ToDOF m_pointDOF[3];
	mutable size_t m_idx_constr;
	/** Pointers to entries in the sparse Jacobian dPhi_dq */
	mutable double *dPhi_dx[3], *dPhi_dy[3];

	/** Pointers to entries in the sparse Jacobian \dot{dPhi_dq} */
	mutable double *dot_dPhi_dx[3], *dot_dPhi_dy[3];

	/** Pointers to entries in the sparse Jacobian d(Phiq*dq)_dq */
	mutable double *dPhiqdq_dx = nullptr, *dPhiqdq_dy = nullptr,
				   *dPhiqdq_dx0 = nullptr, *dPhiqdq_dy0 = nullptr,
				   *dPhiqdq_dx1 = nullptr, *dPhiqdq_dy1 = nullptr;
};

}  // namespace mbse
