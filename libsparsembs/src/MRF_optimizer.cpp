#include <sparsembs/MRF_optimizer.h>

#if defined(SPARSEMBS_HAVE_GTSAM)
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#endif

using namespace sparsembs;
using namespace Eigen;

void sparsembs::MRF_optimizer(
	const MRF_Optimize_Inputs& in, MRF_Optimize_Results& results)
{
#if defined(SPARSEMBS_HAVE_GTSAM)
	MRPT_TRY_START

	gtsam::NonlinearFactorGraph graph;

	gtsam::LevenbergMarquardtOptimizer solver;

	MRPT_TRY_END
#else
	throw std::runtime_error("Requires GTSAM");
#endif
}
