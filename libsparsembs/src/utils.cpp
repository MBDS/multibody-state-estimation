#include <sparsembs/sparsembs-common.h>
#include <sparsembs/sparsembs-utils.h>
#include <mrpt/system/os.h>

using namespace sparsembs;
using namespace Eigen;

mrpt::utils::CTimeLogger sparsembs::timelog;

bool sparsembs::save_matrix(
	cholmod_sparse* tri, const char* filename, cholmod_common* c)
{
	FILE* f = mrpt::system::os::fopen(filename, "wt");
	if (!f) return false;

	const int ret = cholmod_write_sparse(f, tri, NULL, NULL, c);

	fclose(f);

	return ret >= 0;
}

bool sparsembs::save_matrix_dense(
	cholmod_sparse* tri, const char* filename, cholmod_common* c)
{
	FILE* f = mrpt::system::os::fopen(filename, "wt");
	if (!f) return false;

	cholmod_dense* d = cholmod_sparse_to_dense(tri, c);

	const int ret = cholmod_write_dense(f, d, NULL, c);

	cholmod_free_dense(&d, c);

	fclose(f);

	return ret >= 0;
}
