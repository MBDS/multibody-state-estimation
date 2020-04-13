/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/MRF_optimizer.h>

#if defined(SPARSEMBS_HAVE_GTSAM)
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#endif

using namespace mbse;
using namespace Eigen;

void mbse::MRF_optimizer(
	const MRF_Optimize_Inputs& in, MRF_Optimize_Results& results)
{
#if defined(SPARSEMBS_HAVE_GTSAM)
	MRPT_TRY_START

	gtsam::NonlinearFactorGraph graph;
	//	gtsam::LevenbergMarquardtOptimizer solver();

	MRPT_TRY_END
#else
	throw std::runtime_error("Requires GTSAM");
#endif
}
