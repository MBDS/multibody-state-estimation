# ----------------------------------------------------------------------------
#  Auxiliary static library: gtest (Google test framework)
#   Used for Unit Testing
# ----------------------------------------------------------------------------

project(gtest)

include_directories("fused-src")

add_library(gtest 
	STATIC
	fused-src/gtest/gtest-all.cc 
	)
set_target_properties(gtest PROPERTIES FOLDER "Unit tests")

if(CMAKE_COMPILER_IS_GNUCXX AND UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(CMAKE_COMPILER_IS_GNUCXX AND UNIX)

