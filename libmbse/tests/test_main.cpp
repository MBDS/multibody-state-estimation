#include <gtest/gtest.h>

#include <mbse/mbse.h>

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace mbse;

int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);

	return RUN_ALL_TESTS();
}
