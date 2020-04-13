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
#include "mbse-utils.h"
#include "CModelDefinition.h"
#include "CAssembledRigidModel.h"
#include "dynamic-simulators.h"
#include "virtual-sensors.h"

namespace mbse
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

}  // namespace mbse
