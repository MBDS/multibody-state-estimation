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
#include <mbse/Body.h>
#include <mbse/constraints/ConstraintBase.h>
#include <mbse/mbse-utils.h>
#include <mrpt/opengl/CSetOfObjects.h>
#include <functional>
#include <mrpt/core/optional_ref.h>
#include <mrpt/containers/yaml.h>

namespace mbse
{
struct TSymbolicAssembledModel;	 // Frwd. decl.

/** A MBS preprocessed and ready for kinematic/dynamic simulations. */
class AssembledRigidModel;

/** The class for user-defined MBS problems.
 *  After construction, the user must ensure that all fields are correctly
 * filled-in to define the MBS joints & elements.
 * Relative coordinates can be directly added into rDOFs_.
 *
 * Mechanism models can be more easily defined in a YAML file and loaded via
 * FromYAML().
 *
 */
class ModelDefinition
{
   public:
	/** Constructor: creates an empty model */
	ModelDefinition() = default;

	/** Builds a multibody problem from a YAML definition, in the format shown
	 * in examples. */
	static ModelDefinition FromYAML(const mrpt::containers::yaml& c);

	/** Completely erases all defined points, joints, bodies, parameters, etc of
	 * this object and leaves it blank. */
	void clear();

	/** Set the number of points (all, fixed and variables) in the MBS problem
	 */
	void setPointCount(const size_t n) { points_.resize(n); }

	/** Get the number of points (all, fixed and variables) in the MBS problem
	 */
	size_t getPointCount() const { return points_.size(); }

	/** Set the (0-based) i'th point's coordinates and its type (fixed/variable)
	 *  TODO: Initial position problem will refine these positions if needed.
	 */
	void setPointCoords(
		const size_t i, const mrpt::math::TPoint2D& coords,
		const bool is_fixed = false);

	const Point2& getPointInfo(const size_t i) const
	{
		ASSERT_(i < points_.size());
		return points_[i];
	}

	/** Appends a new empty body to the end of bodies_ and returns a reference
	 * for further fill-in. If no name is provided, an automatic numbered name
	 * will be generated. The user must fill in the returned object with the
	 * desired values.
	 */
	Body& addBody(const std::string& name = std::string(""));

	/** Introduces a new constraint in the MBS.
	 *  Note that fixed-distance constraints due to rigid bodies don't have to
	 * be added by the user, since they're always automatically added. Check
	 * derived classes of ConstraintBase to see the list of all possible
	 * constraints.
	 * \note Expected arguments are those of the constraint class ctor.
	 */
	template <class CONSTRAINT_CLASS, typename... _Args>
	void addConstraint(_Args&&... args)
	{
		// Make a copy as dynamic memory and save as a smart pointer
		// of the base class.
		constraints_.emplace_back(std::make_shared<CONSTRAINT_CLASS>(args...));
	}

	/** Process the MBS definitions and assemble all the required symbolic
	 * structures to enable kinematic/dynamic simulations of the MBS.
	 * \param[out] out_armi Must be created with *this as parent model.
	 *  Returns the "list of instructions" to build an actual assembled model of
	 * class AssembledRigidModel.
	 */
	void assembleRigidMBS(TSymbolicAssembledModel& out_armi) const;

	/** Process the MBS definitions and assemble all the required symbolic
	 * structures to enable kinematic/dynamic simulations of the MBS. Invoke
	 * methods of the returned object.
	 * Optionally, relative coordinates may be added.
	 */
	std::shared_ptr<AssembledRigidModel> assembleRigidMBS() const;

	/** Additional relative coordinates */
	std::vector<RelativeDOF> rDOFs_;

	const std::vector<ConstraintBase::Ptr>& constraints() const
	{
		return constraints_;
	}
	const std::vector<Body>& bodies() const { return bodies_; }

	std::vector<ConstraintBase::Ptr>& constraints() { return constraints_; }
	std::vector<Body>& bodies() { return bodies_; }

   protected:
	/** @name Data
		@{ */
	std::vector<Point2> points_;  //!< ALL points (fixed and variables)
	std::vector<Body> bodies_;	//!< Bodies

	/** The list of all constraints (of different kinds/classes).
	 * \note Constant-distance constraints for rigid bodies are NOT included in
	 * this list, but are automatically added to constraints list in ARM.
	 */
	std::vector<ConstraintBase::Ptr> constraints_;

	/** @} */  // end data --------------

	mutable bool already_added_fixed_len_constraints_ = false;

};	// end class ModelDefinition

}  // namespace mbse
