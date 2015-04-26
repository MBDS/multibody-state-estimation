#include <gtest/gtest.h>

#include <sparsembs/sparsembs.h>

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace sparsembs;

int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
