# ----------------------------------------------------------------------------
#   TESTS
# ----------------------------------------------------------------------------
PROJECT(tests)

# Fix a "bug" in VS11 (MSVC 2012):
if(MSVC)
	add_definitions(-D_VARIADIC_MAX=10)
endif(MSVC)


# Compile gtest from sources:
add_subdirectory(gtest-1.6.0)
INCLUDE_DIRECTORIES("gtest-1.6.0/fused-src")


SET(EXECUTABLE_OUTPUT_PATH "${${PROJECT_NAME}_BINARY_DIR}")

# -----------------------------
# Add the "make test" target:
# -----------------------------
ADD_CUSTOM_TARGET(test)
set_target_properties(test PROPERTIES FOLDER "Unit tests")


# Tests based on Google gtest:
# -----------------------------

# If using system library, add C++ flags:
IF (UNIX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
ENDIF (UNIX)

FILE(GLOB TEST_SRCS "*.cpp")

# Test project:
ADD_EXECUTABLE(unit_tests ${TEST_SRCS})
TARGET_LINK_LIBRARIES(unit_tests ${MRPT_LIBS} ${SPARSEMBS_LIBS})
IF(NOT "${DEPENCENCIES_SPARSEMBS}" STREQUAL "")
	ADD_DEPENDENCIES(unit_tests ${DEPENCENCIES_SPARSEMBS})
ENDIF(NOT "${DEPENCENCIES_SPARSEMBS}" STREQUAL "")

# Add the required libraries for linking:
TARGET_LINK_LIBRARIES(unit_tests gtest)
set_target_properties(unit_tests PROPERTIES FOLDER "Unit tests")

# Run it:
GET_TARGET_PROPERTY(GENERATED_EXE unit_tests LOCATION)
ADD_CUSTOM_TARGET(run_tests COMMAND "${GENERATED_EXE}" "${${PROJECT_NAME}__SOURCE_DIR}")
set_target_properties(run_tests PROPERTIES FOLDER "Unit tests")
ADD_DEPENDENCIES(run_tests unit_tests)
ADD_DEPENDENCIES(test run_tests)
