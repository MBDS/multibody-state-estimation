/*+-------------------------------------------------------------------------+
  |            Multi Body State Estimation (mbse) C++ library               |
  |                                                                         |
  | Copyright (C) 2014-2020 University of Almeria                           |
  | Copyright (C) 2020 University of Salento                                |
  | See README for list of authors and papers                               |
  | Distributed under 3-clause BSD license                                  |
  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
  +-------------------------------------------------------------------------+ */

#include <mbse/mbse-common.h>
#include <mbse/mbse-utils.h>
#include <mrpt/system/os.h>

using namespace mbse;
using namespace Eigen;

CTimeLogger mbse::timelog;

bool mbse::save_matrix(
	cholmod_sparse* tri, const char* filename, cholmod_common* c)
{
	FILE* f = mrpt::system::os::fopen(filename, "wt");
	if (!f) return false;

	const int ret = cholmod_write_sparse(f, tri, NULL, NULL, c);

	fclose(f);

	return ret >= 0;
}

bool mbse::save_matrix_dense(
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
