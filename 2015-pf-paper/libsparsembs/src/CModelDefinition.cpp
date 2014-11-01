#include <sparsembs/CModelDefinition.h>
#include <sparsembs/CAssembledModelRigid.h>
#include <mrpt/opengl.h>


using namespace sparsembs;
using namespace std;

CModelDefinition::CModelDefinition() :
	m_already_added_fixed_len_constraints(false)
{
}


CBody & CModelDefinition::addBody(const std::string &name )
{
	ASSERTMSG_(!m_already_added_fixed_len_constraints,"Can't modify model after assembling!")

	// Build name if none provided:
	const std::string nam = name.empty() ?
		mrpt::format("body%u",static_cast<unsigned int>(m_bodies.size()) )
		:
		name;

	// Create:
	m_bodies.resize(m_bodies.size() +1);

	// Set name & return:
	CBody & new_body = *m_bodies.rbegin();
	new_body.name = nam;

	return new_body;
}


/** Completely erases all defined points, joints, bodies, parameters, etc of this object and leaves it blank. */
void CModelDefinition::clear()
{
	*this = CModelDefinition();
}

MRPT_TODO("Initial position problem should refine these positions if needed.")

void CModelDefinition::setPointCoords(const size_t i, const TPoint2D &coords, const bool is_fixed)
{
	ASSERT_(i<m_points.size())

	TMBSPoint & pt = m_points[i];
	pt.coords = coords;
	pt.fixed = is_fixed;
}


// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
void CModelDefinition::assembleRigidMBS(TSymbolicAssembledModel &armi) const
{
	armi.clear();

	// 1) Count number of natural coordinates which are unknowns ==> # of scalar unknowns in vector q
	// ---------------------------------------------------------------------------------------------------
	const size_t nPts = m_points.size();
	for (size_t i=0;i<nPts;i++)
	{
		const TMBSPoint &pt = m_points[i];
		// If it's not fixed, add to the list of DOFs
		if (!pt.fixed)
		{
			armi.DOFs.push_back( TDOF(i,0) ); // x
			armi.DOFs.push_back( TDOF(i,1) ); // y
		}
	}

	// 2) For each rigid body, automatically add one constant-distance constraint:
	// ---------------------------------------------------------------------------------------------------
	if (!m_already_added_fixed_len_constraints)
	{
		for (size_t i=0;i<m_bodies.size();i++)
		{
			const CBody &b = m_bodies[i];
			ASSERT_(b.points[0]<m_points.size())
			ASSERT_(b.points[1]<m_points.size())

			const_cast<CModelDefinition*>(this)->addConstraint( CConstraintConstantDistance( b.points[0],b.points[1],b.length() ) );
		}
		// Mark these constraints as added:
		m_already_added_fixed_len_constraints = true;
	}
}



// -------------------------------------------------------------------
//                      assembleRigidMBS
// -------------------------------------------------------------------
CAssembledRigidModelPtr CModelDefinition::assembleRigidMBS()
{
	// 1) Build "symbolic" assembly:
	TSymbolicAssembledModel armi(*this);
	this->assembleRigidMBS(armi);
	
	// 2) Actual assembly:
	return CAssembledRigidModelPtr(new CAssembledRigidModel(armi));
}

