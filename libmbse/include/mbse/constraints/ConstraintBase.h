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

#include <memory>
#include <mrpt/opengl/CRenderizable.h>

namespace mbse
{
class AssembledRigidModel;

/** The virtual base class of all constraint types. */
class ConstraintBase
{
   public:
	/** A smart pointer type for constraints */
	using Ptr = std::shared_ptr<ConstraintBase>;

	/** Alloc space for the needed rows, columns, etc. in the constraint vector
	 * and its Jacobians Each class should save a reference/pointer to the newly
	 * created elements, so they can be quickly updated in the \a update()
	 * method when invoked in the future. This method is called after the MBS is
	 * completely defined by the user, and before starting any kinematic/dynamic
	 * simulation.
	 */
	virtual void buildSparseStructures(AssembledRigidModel& arm) const = 0;

	/** Update the previously allocated values given the current state of the
	 * MBS. This is called a very large number of times during simulations. */
	virtual void update(AssembledRigidModel& arm) const = 0;

	/** Virtual destructor (required in any virtual base) */
	virtual ~ConstraintBase();

	/** Clone operator for smart pointers */
	virtual Ptr clone() const = 0;

	/** Creates a 3D representation of the constraint, if applicable (e.g.the
	 * line of a fixed slider).
	 * \return false if the constraint has no 3D representation
	 */
	virtual bool get3DRepresentation(
		[[maybe_unused]]  //
		mrpt::opengl::CRenderizable::Ptr& inout_obj) const
	{
		return false;
	}
};
}  // namespace mbse
