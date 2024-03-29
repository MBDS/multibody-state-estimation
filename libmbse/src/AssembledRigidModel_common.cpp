/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2024 University of Almeria                           |
  | Copyright (C) 2021 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/AssembledRigidModel.h>
#include <mbse/constraints/ConstraintRelativeAngle.h>
#include <mbse/constraints/ConstraintRelativeAngleAbsolute.h>
#include <mrpt/opengl.h>
#include <iostream>

using namespace mbse;
using namespace Eigen;
using namespace mrpt::math;
using namespace mrpt;
using namespace std;

const double DEFAULT_GRAVITY[3] = {0, -9.81, 0};

/** Constructor */
AssembledRigidModel::AssembledRigidModel(const TSymbolicAssembledModel& armi)
	: mechanism_(armi.model)
{
	for (int i = 0; i < 3; i++) gravity_[i] = DEFAULT_GRAVITY[i];

	const auto nEuclideanDOFs = armi.DOFs.size();
	const auto nRelativeDOFs = armi.rDOFs.size();

	const auto nDOFs = nEuclideanDOFs + nRelativeDOFs;

	ASSERTMSG_(
		nEuclideanDOFs > 0,
		"Trying to assemble model with 0 Natural Coordinate DOFs");

	q_.setZero(nDOFs);
	dotq_.setZero(nDOFs);
	ddotq_.setZero(nDOFs);
	Q_.setZero(nDOFs);

	// Keep a copy of DOF info:
	DOFs_ = armi.DOFs;
	rDOFs_ = armi.rDOFs;

	// Build reverse-lookup table & initialize initial values for "q":
	points2DOFs_.resize(armi.model.getPointCount());

	for (dof_index_t i = 0; i < nEuclideanDOFs; i++)
	{
		const size_t pt_idx = DOFs_.at(i).point_index;
		const Point2& pt = mechanism_.getPointInfo(pt_idx);
		switch (DOFs_.at(i).point_dof)
		{
			case PointDOF::X:
				q_[i] = pt.coords.x;
				points2DOFs_.at(pt_idx).dof_x = i;
				break;
			case PointDOF::Y:
				q_[i] = pt.coords.y;
				points2DOFs_.at(pt_idx).dof_y = i;
				break;
			// case 2: //z
			//	q_[i] = pt.coords.z;
			//	break;
			default:
				THROW_EXCEPTION("Unexpected value for DOFs_[i].point_dof");
		};
	}

	// Save the number of DOFs as the number of columsn in sparse Jacobians:
	defineSparseMatricesColumnCount(nDOFs);

	// Generate constraint equations & create structure of sparse Jacobians:
	// ----------------------------------------------------------------------
	const std::vector<ConstraintBase::Ptr>& parent_constraints =
		mechanism_.constraints();

	// 1/2: Constraints
	const size_t nConst = parent_constraints.size();
	constraints_.resize(nConst);
	for (size_t i = 0; i < nConst; i++)
	{
		constraints_[i] = parent_constraints[i]->clone();
	}

	// 2/2: Constraints from relative coordinates:
	relCoordinate2Index_.resize(rDOFs_.size());
	for (size_t i = 0; i < rDOFs_.size(); i++)
	{
		const auto& relConstr = rDOFs_[i];
		const dof_index_t idxInQ = nEuclideanDOFs + i;

		if (std::holds_alternative<RelativeAngleDOF>(relConstr))
		{
			const auto& c = std::get<RelativeAngleDOF>(relConstr);

			// Add constraint:
			auto co = std::make_shared<ConstraintRelativeAngle>(
				c.point_idx0, c.point_idx1, c.point_idx2, idxInQ);
			constraints_.push_back(co);

			// reverse look up table:
			relCoordinate2Index_[i] = idxInQ;
		}
		else if (std::holds_alternative<RelativeAngleAbsoluteDOF>(relConstr))
		{
			const auto& c = std::get<RelativeAngleAbsoluteDOF>(relConstr);

			// Add constraint:
			auto co = std::make_shared<ConstraintRelativeAngleAbsolute>(
				c.point_idx0, c.point_idx1, idxInQ);
			constraints_.push_back(co);

			// reverse look up table:
			relCoordinate2Index_[i] = idxInQ;

			// Initial value from coordinates:
			const auto pt0 = getPointCurrentCoords(c.point_idx0);
			const auto pt1 = getPointCurrentCoords(c.point_idx1);
			const auto ptV = pt1 - pt0;
			if (ptV.norm() != 0)
			{
				auto expectedAng = std::atan2(ptV.y, ptV.x);
				if (std::abs(q_[idxInQ] - expectedAng) > 0.01)
				{
					printf(
						"*WARNING* RelativeAngleAbsoluteDOF at q[%i]=%f deg: "
						"Overriding with angle from coordinates = %f deg.\n",
						static_cast<int>(idxInQ), mrpt::RAD2DEG(q_[idxInQ]),
						mrpt::RAD2DEG(expectedAng));
					q_[idxInQ] = expectedAng;
				}
			}
		}
		else
		{
			THROW_EXCEPTION("Unknown type of relative coordinate");
		}
	}

	// Final step: build structures
	for (auto& c : constraints_) c->buildSparseStructures(*this);
}

void AssembledRigidModel::getGravityVector(
	double& gx, double& gy, double& gz) const
{
	gx = gravity_[0];
	gy = gravity_[1];
	gz = gravity_[2];
}

void AssembledRigidModel::setGravityVector(
	const double gx, const double gy, const double gz)
{
	gravity_[0] = gx;
	gravity_[1] = gy;
	gravity_[2] = gz;
}

/** Call all constraint objects and command them to update their corresponding
 * parts in the sparse Jacobians */
void AssembledRigidModel::update_numeric_Phi_and_Jacobians()
{
	// Update numeric values of the constraint Jacobians:
	for (size_t i = 0; i < constraints_.size(); i++)
		constraints_[i]->update(*this);
}

void AssembledRigidModel::realize_operating_point() const
{
	// Update numeric values of the constraint Jacobians:
	for (size_t i = 0; i < constraints_.size(); i++)
		constraints_[i]->realizeOperatingPoint(*this);
}

/** Returns a 3D visualization of the model */
void AssembledRigidModel::getAs3DRepresentation(
	mrpt::opengl::CSetOfObjects::Ptr& outObj,
	const Body::TRenderParams& rp) const
{
	if (!outObj)
		outObj = mrpt::opengl::CSetOfObjects::Create();
	else
		outObj->clear();

	// Render constraints:
	const std::vector<ConstraintBase::Ptr>& parent_constr =
		mechanism_.constraints();

	const size_t nConstr = parent_constr.size();
	for (size_t i = 0; i < nConstr; i++)
	{
		const ConstraintBase* constr_ptr = parent_constr[i].get();
		mrpt::opengl::CRenderizable::Ptr gl_obj;
		if (constr_ptr->get3DRepresentation(gl_obj))
		{
			// Insert in 3D scene:
			outObj->insert(gl_obj);
		}
	}

	// Render "ground" points:
	if (rp.show_grounds)
	{
		const size_t nPts = mechanism_.getPointCount();
		for (size_t i = 0; i < nPts; i++)
		{
			const Point2& pt = mechanism_.getPointInfo(i);
			if (!pt.fixed) continue;
			// This is a fixed point:
			outObj->insert(this->internal_render_ground_point(pt, rp));
		}
	}

	// Render bodies:
	const std::vector<Body>& parent_bodies = mechanism_.bodies();
	const size_t nBodies = parent_bodies.size();
	gl_objects_.resize(nBodies);

	for (size_t i = 0; i < nBodies; i++)
	{
		const Body& b = parent_bodies[i];

		mrpt::opengl::CRenderizable::Ptr gl_obj = b.get3DRepresentation();

		// Insert in 3D scene:
		outObj->insert(gl_obj);

		// And save reference for quickly update the 3D pose in the future
		// during animations:
		gl_objects_[i] = gl_obj;
	}

	// Place each body in its current pose:
	this->update3DRepresentation(rp);
}

/** Animates a 3D representation of the MBS, previously built in
 * getAs3DRepresentation() \sa getAs3DRepresentation
 */
void AssembledRigidModel::update3DRepresentation(
	const Body::TRenderParams& rp) const
{
	const std::vector<Body>& parent_bodies = mechanism_.bodies();

	const size_t nBodies = parent_bodies.size();
	if (gl_objects_.size() != nBodies)
	{
		std::cerr << "[AssembledRigidModel::update3DRepresentation] Warning: "
					 "Opengl model is not initialized.\n";
		return;
	}

	for (size_t i = 0; i < nBodies; i++)
	{
		mrpt::opengl::CRenderizable::Ptr& obj = gl_objects_[i];
		ASSERT_(obj);

		// Recover the 2D pose from 2 points:
		const Body& b = parent_bodies[i];
		const size_t i0 = b.points[0];
		const size_t i1 = b.points[1];

		const Point2& pnt0 = mechanism_.getPointInfo(i0);
		const Point2& pnt1 = mechanism_.getPointInfo(i1);

		const Point2ToDOF dof0 = points2DOFs_[i0];
		const Point2ToDOF dof1 = points2DOFs_[i1];

		const double& p0x =
			(dof0.dof_x != INVALID_DOF) ? q_[dof0.dof_x] : pnt0.coords.x;
		const double& p0y =
			(dof0.dof_y != INVALID_DOF) ? q_[dof0.dof_y] : pnt0.coords.y;

		const double& p1x =
			(dof1.dof_x != INVALID_DOF) ? q_[dof1.dof_x] : pnt1.coords.x;
		const double& p1y =
			(dof1.dof_y != INVALID_DOF) ? q_[dof1.dof_y] : pnt1.coords.y;

		const double theta = atan2(p1y - p0y, p1x - p0x);

		obj->setPose(
			mrpt::poses::CPose3D(p0x, p0y, 0, theta, DEG2RAD(0), DEG2RAD(0)));

		// Update transparency:
		if (rp.render_style == Body::reLine)
		{
			auto gl_line =
				mrpt::ptr_cast<mrpt::opengl::CSetOfObjects>::from(obj)
					->getByClass<mrpt::opengl::CSimpleLine>();
			if (gl_line) gl_line->setColorA_u8(rp.line_alpha);
		}
	}
}

size_t AssembledRigidModel::addNewRowToConstraints()
{
	const size_t idx = Phi_.size();
	const size_t m = idx + 1;  // new size
	resizeConstraintCount(m);
	return idx;
}

/** Only to be called between objects created from the same symbolic model, this
 * method replicates the state of "o" into "this". */
void AssembledRigidModel::copyStateFrom(const AssembledRigidModel& o)
{
#ifdef _DEBUG
	// Security checks, just in case...
	ASSERT_(this->q_.size() == o.q_.size());
	ASSERT_(this->dotq_.size() == o.dotq_.size());
	ASSERT_(this->DOFs_.size() == o.DOFs_.size());
	ASSERT_(this->Phi_.size() == o.Phi_.size());
	const double* ptr_q0 = &q_[0];
	const double* ptr_dotq0 = &dotq_[0];
#endif

	this->q_ = o.q_;
	this->dotq_ = o.dotq_;

#ifdef _DEBUG
	ASSERT_(
		ptr_q0 == &q_[0]);	// make sure the vectors didn't suffer mem
							// reallocation, since we save pointers to these!
	ASSERT_(
		ptr_dotq0 ==
		&dotq_[0]);	 // make sure the vectors didn't suffer mem reallocation,
					 // since we save pointers to these!
#endif
}

/** Copies the opengl object from another instance */
void AssembledRigidModel::copyOpenGLRepresentationFrom(
	const AssembledRigidModel& o)
{
	this->gl_objects_ = o.gl_objects_;
}

/** Retrieves the current coordinates of a point, which may include either fixed
 * or variable components */
void AssembledRigidModel::getPointCurrentCoords(
	const size_t pt_idx, mrpt::math::TPoint2D& pt) const
{
	const Point2& pt_info = mechanism_.getPointInfo(pt_idx);
	const Point2ToDOF& pt_dofs = points2DOFs_[pt_idx];

	pt.x =
		(pt_dofs.dof_x != INVALID_DOF) ? q_[pt_dofs.dof_x] : pt_info.coords.x;
	pt.y =
		(pt_dofs.dof_y != INVALID_DOF) ? q_[pt_dofs.dof_y] : pt_info.coords.y;
}

/** Retrieves the current velocity of a point, which may include either fixed or
 * variable components */
void AssembledRigidModel::getPointCurrentVelocity(
	const size_t pt_idx, mrpt::math::TPoint2D& vel) const
{
	const Point2ToDOF& pt_dofs = points2DOFs_[pt_idx];

	vel.x = (pt_dofs.dof_x != INVALID_DOF) ? dotq_[pt_dofs.dof_x] : 0;
	vel.y = (pt_dofs.dof_y != INVALID_DOF) ? dotq_[pt_dofs.dof_y] : 0;
}

/** Computes the current coordinates of a point fixed to a given body, given its
 * relative coordinates wrt to system X:pt0->pt1, Y: orthogonal */
void AssembledRigidModel::getPointOnBodyCurrentCoords(
	const size_t body_index, const mrpt::math::TPoint2D& relative_pt,
	mrpt::math::TPoint2D& out_pt) const
{
	ASSERTDEB_(body_index < mechanism_.bodies().size());

	const Body& b = mechanism_.bodies()[body_index];

	mrpt::math::TPoint2D q[2];
	this->getPointCurrentCoords(b.points[0], q[0]);
	this->getPointCurrentCoords(b.points[1], q[1]);

	const double L = b.length();
	ASSERTDEB_(L > 0);
	const double Linv = 1.0 / L;

	mrpt::math::TPoint2D u, v;	// unit vectors in X,Y,Z local to the body

	u = (q[1] - q[0]) * Linv;
	v.x = -u.y;
	v.y = u.x;

	out_pt.x = q[0].x + u.x * relative_pt.x + v.x * relative_pt.y;
	out_pt.y = q[0].y + u.y * relative_pt.x + v.y * relative_pt.y;
}

/* Render a ground point */
mrpt::opengl::CSetOfObjects::Ptr
	AssembledRigidModel::internal_render_ground_point(
		const Point2& pt, const Body::TRenderParams& rp) const
{
	mrpt::opengl::CSetOfObjects::Ptr obj =
		mrpt::opengl::CSetOfObjects::Create();
	obj->setLocation(pt.coords.x, pt.coords.y, 0);

	const double support_LXZ = 0.03;
	const double support_LY = 0.05;

	auto gl_box = mrpt::opengl::CBox::Create(
		mrpt::math::TPoint3D(
			-0.5 * support_LXZ, -support_LY, -0.5 * support_LXZ),
		mrpt::math::TPoint3D(0.5 * support_LXZ, 0, 0.5 * support_LXZ), false);
	gl_box->setColor(0, 0, 0.7);
	obj->insert(gl_box);

	return obj;
}

/** Evaluate current energy of the system. */
void AssembledRigidModel::evaluateEnergy(
	AssembledRigidModel::TEnergyValues& e) const
{
	timelog().enter("evaluateEnergy");

	e = TEnergyValues();  // Reset to zero

	const std::vector<Body>& bodies = mechanism_.bodies();
	for (size_t i = 0; i < bodies.size(); i++)
	{
		const Body& b = bodies[i];

		mrpt::math::TPoint2D dq[2];
		this->getPointCurrentVelocity(b.points[0], dq[0]);
		this->getPointCurrentVelocity(b.points[1], dq[1]);

		const Matrix<double, 2, 1> dq0 = Matrix<double, 2, 1>(dq[0].x, dq[0].y);
		const Matrix<double, 2, 1> dq1 = Matrix<double, 2, 1>(dq[1].x, dq[1].y);

		const Matrix2d& M00 = b.getM00();
		const Matrix2d& M01 = b.getM01();
		const Matrix2d& M11 = b.getM11();

		// Kinetic energy: 0.5 * [q0 q1] * [M00 M10;M10' M11] * [q0 q1]'
		e.E_kin +=
			0.5 * (dq0.transpose() * M00 * dq0 + dq1.transpose() * M11 * dq1)
					  .coeff(0, 0) +
			(dq0.transpose() * M01 * dq1).coeff(0, 0);

		// Potential energy:
		mrpt::math::TPoint2D
			global_cog;	 // current COG position, in global coords:
		this->getPointOnBodyCurrentCoords(i, b.cog(), global_cog);

		e.E_pot -= b.mass() * (this->gravity_[0] * global_cog.x +
							   this->gravity_[1] * global_cog.y +
							   this->gravity_[2] * 0 /*global_cog.z*/);
	}

	e.E_total = e.E_kin + e.E_pot;

	timelog().leave("evaluateEnergy");
}

void AssembledRigidModel::printCoordinates(std::ostream& o) const
{
	MRPT_START

	ASSERT_EQUAL_(static_cast<size_t>(q_.size()), DOFs_.size() + rDOFs_.size());

	o << "[AssembledRigidModel] |q|=" << q_.size() << ", " << DOFs_.size()
	  << " natural, " << rDOFs_.size() << " relative coordinates.\n";

	o << "Natural coordinates:\n";
	for (size_t i = 0; i < DOFs_.size(); i++)
	{
		o << " q[" << i << "]: " << dof2letter(DOFs_[i].point_dof)
		  << DOFs_[i].point_index << "\n";
	}
	if (!rDOFs_.empty())
	{
		o << "Relative coordinates:\n";
		for (size_t i = 0; i < rDOFs_.size(); i++)
		{
			o << " q[" << (i + DOFs_.size()) << "]: ";
			const auto& relConstr = rDOFs_[i];
			if (std::holds_alternative<RelativeAngleAbsoluteDOF>(relConstr))
			{
				const auto& c = std::get<RelativeAngleAbsoluteDOF>(relConstr);
				o << "relativeAngleWrtGround(" << c.point_idx0 << " - "
				  << c.point_idx1 << ")";
			}
			else
			{
				o << "???";
			}
			o << "\n";
		}
	}

	MRPT_END
}

void AssembledRigidModel::printConstraints(std::ostream& o) const
{
	MRPT_START

	o << "[AssembledRigidModel] m=" << constraints_.size() << " constraints.\n";
	for (size_t i = 0; i < constraints_.size(); i++)
	{
		o << "- constraint[" << i << "]:\n";
		constraints_.at(i)->print(o);
	}

	MRPT_END
}
