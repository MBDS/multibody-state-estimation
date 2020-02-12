#pragma once

#include "sparsembs-common.h"
#include "sparsembs-utils.h"
#include "CModelDefinition.h"
#include "CAssembledRigidModel.h"
#include "dynamic-simulators.h"
#include "virtual-sensors.h"

namespace sparsembs
{
using namespace Eigen;
using namespace mrpt::math;

struct MRF_Optimize_Inputs
{
	MRF_Optimize_Inputs(CDynamicSimulatorBase& d) : dysim(d) {}

	CDynamicSimulatorBase& dysim;
};

struct MRF_Optimize_Results
{
};

/** Runs a batch optimization process
 */
void MRF_optimizer(
	const MRF_Optimize_Inputs& in, MRF_Optimize_Results& results);

}  // namespace sparsembs
