#pragma once

#include "sparsembs-common.h"
#include "CBody.h"
#include "constraints.h"
#include <functional>
#include "sparsembs-utils.h"
#include <Eigen/UmfPackSupport>
#include <mrpt/opengl/CSetOfObjects.h>

namespace sparsembs
{
using namespace Eigen;
using namespace mrpt::math;

struct TSymbolicAssembledModel;  // Frwd. decl.

class CAssembledRigidModel;  //!< A MBS preprocessed and ready for
							 //!< kinematic/dynamic simulations.
typedef stlplus::smart_ptr_clone<CAssembledRigidModel>
	CAssembledRigidModelPtr;  //!< A smart-pointer for CAssembledRigidModel

/** The class for user-defined MBS problems.
 *  After construction, the user must assure that all fields are correctly
 * filled-in to define the MBS joints & elements.
 */
class CModelDefinition
{
   public:
	/** Constructor: creates an empty model */
	CModelDefinition();

	/** Completely erases all defined points, joints, bodies, parameters, etc of
	 * this object and leaves it blank. */
	void clear();

	/** Set the number of points (all, fixed and variables) in the MBS problem
	 */
	void setPointCount(const size_t n) { m_points.resize(n); }

	/** Get the number of points (all, fixed and variables) in the MBS problem
	 */
	size_t getPointCount() const { return m_points.size(); }

	/** Set the (0-based) i'th point's coordinates and its type (fixed/variable)
	 *  TODO: Initial position problem will refine these positions if needed.
	 */
	void setPointCoords(
		const size_t i, const TPoint2D& coords, const bool is_fixed = false);

	const TMBSPoint& getPointInfo(const size_t i) const
	{
		ASSERTDEB_(i < m_points.size()) return m_points[i];
	}

	/** Appends a new empty body to the end of m_bodies and returns a reference
	 * for further fill-in. If no name is provided, an automatic numbered name
	 * will be generated. The user must fill in the returned object with the
	 * desired values.
	 */
	CBody& addBody(const std::string& name = std::string(""));

	/** Introduces a new constraint in the MBS.
	 *  Note that fixed-distance constraints due to rigid bodies don't have to
	 * be added by the user, since they're always automatically added. Check
	 * derived classes of CConstraintBase to see the list of all possible
	 * constraints. \note A copy is made from the passed object, so it can be
	 * safely deleted upon return.
	 */
	template <class CONSTRAINT_CLASS>
	void addConstraint(const CONSTRAINT_CLASS& c)
	{
		m_constraints.push_back(CConstraintBasePtr(new CONSTRAINT_CLASS(
			c)));  // Make a copy as dynamic memory and save as a smart pointer
				   // of the base class.
	}

	/** Process the MBS definitions and assemble all the required symbolic
	 * structures to enable kinematic/dynamic simulations of the MBS.
	 * \param[out] out_armi Must be created with *this as parent model.
	 *  Returns the "list of instructions" to build an actual assembled model of
	 * class CAssembledRigidModel.
	 */
	void assembleRigidMBS(TSymbolicAssembledModel& out_armi) const;

	/** Process the MBS definitions and assemble all the required symbolic
	 * structures to enable kinematic/dynamic simulations of the MBS. Invoke
	 * methods of the returned object.
	 */
	CAssembledRigidModelPtr assembleRigidMBS();

	const std::vector<CConstraintBasePtr>& getConstraints() const
	{
		return m_constraints;
	}
	const std::vector<CBody>& getBodies() const { return m_bodies; }

	std::vector<CConstraintBasePtr>& getConstraints() { return m_constraints; }
	std::vector<CBody>& getBodies() { return m_bodies; }

   protected:
	/** @name Data
		@{ */
	std::vector<TMBSPoint> m_points;  //!< ALL points (fixed and variables)
	std::vector<CBody> m_bodies;  //!< Bodies

	/** The list of all constraints (of different kinds/classes).
	 * \note Constant-distance constraints for rigid bodies are NOT included in
	 * this list, but are automatically added to constraints list in ARM.
	 */
	std::vector<CConstraintBasePtr> m_constraints;

	/** @} */  // end data --------------

	mutable bool m_already_added_fixed_len_constraints;  //!< Self-explanatory

};  // end class CModelDefinition

}  // namespace sparsembs
