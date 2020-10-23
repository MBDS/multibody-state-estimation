# ----------------------------------------------------------------------------
#   TESTS
# ----------------------------------------------------------------------------
project(tests)

# Fix a "bug" in VS11 (MSVC 2012):
if(MSVC)
	add_definitions(-D_VARIADIC_MAX=10)
endif()

# Tests based on Google gtest:
# -----------------------------
# Compile gtest from sources:
add_subdirectory(gtest-1.6.0)
include_directories("gtest-1.6.0/fused-src")

add_custom_target(test_legacy)
add_custom_target(tests_build_all)
set_target_properties(test_legacy PROPERTIES FOLDER "unit tests")
set_target_properties(tests_build_all PROPERTIES FOLDER "unit tests")


# If using system library, add C++ flags:
if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif ()

# Define a test:
function (mbse_define_test NAME)
  add_executable(test_${NAME} test_main.cpp ${NAME}.cpp)

  # Add the required libraries for linking:
  target_link_libraries(test_${NAME}
	gtest
	${MRPT_LIBRARIES}
	mbse::mbse
	)

  set_target_properties(test_${NAME} PROPERTIES FOLDER "Unit tests")

  add_test(test_${NAME}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target test_${NAME})
  add_test(test_${NAME}_run ${EXECUTABLE_OUTPUT_PATH}/test_${NAME})
  set_tests_properties(test_${NAME}_run PROPERTIES DEPENDS test_${NAME}_build)

  # Run it:
  set(GENERATED_EXE "$<TARGET_FILE:test_${NAME}>")
  add_custom_target(run_tests_${NAME} COMMAND "${GENERATED_EXE}")
  add_dependencies(run_tests_${NAME} test_${NAME})
  add_dependencies(test_legacy run_tests_${NAME})
  add_dependencies(tests_build_all test_${NAME})
endfunction()


# List of tests:
mbse_define_test(fourbars)
mbse_define_test(dynamics-solvers)

mbse_define_test(factor-euler-integrator)
mbse_define_test(factor-trapezoidal-integrator)
mbse_define_test(factor-icoords-dynamics-jacobian)
mbse_define_test(factor-dynamics-jacobian)
mbse_define_test(factor-gyroscope-jacobian)

