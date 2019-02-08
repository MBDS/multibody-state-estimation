# ----------------------------------------------------------------------------
#   TESTS
# ----------------------------------------------------------------------------
project(tests)

# Fix a "bug" in VS11 (MSVC 2012):
if(MSVC)
	add_definitions(-D_VARIADIC_MAX=10)
endif(MSVC)


# Compile gtest from sources:
add_subdirectory(gtest-1.6.0)
include_directories("gtest-1.6.0/fused-src")


set(EXECUTABLE_OUTPUT_PATH "${${PROJECT_NAME}_BINARY_DIR}")

# -----------------------------
# Add the "make test" target:
# -----------------------------
add_custom_target(test)
set_target_properties(test PROPERTIES FOLDER "Unit tests")


# Tests based on Google gtest:
# -----------------------------

# If using system library, add C++ flags:
if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif (UNIX)

file(GLOB TEST_SRCS "*.cpp")

# Test project:
add_executable(unit_tests ${TEST_SRCS})
target_link_libraries(unit_tests ${MRPT_LIBS} ${SPARSEMBS_LIBS})
if(NOT "${DEPENCENCIES_SPARSEMBS}" STREQUAL "")
	add_dependencies(unit_tests ${DEPENCENCIES_SPARSEMBS})
endif(NOT "${DEPENCENCIES_SPARSEMBS}" STREQUAL "")

# Add the required libraries for linking:
target_link_libraries(unit_tests gtest)
set_target_properties(unit_tests PROPERTIES FOLDER "Unit tests")

# Run it:
get_target_property(GENERATED_EXE unit_tests LOCATION)
add_custom_target(run_tests COMMAND "${GENERATED_EXE}" "${${PROJECT_NAME}__SOURCE_DIR}")
set_target_properties(run_tests PROPERTIES FOLDER "Unit tests")
add_dependencies(run_tests unit_tests)
add_dependencies(test run_tests)
