project(MBSKIT)

cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)  # Required by CMake 2.7+
endif(COMMAND cmake_policy)

set(LIBRARY_OUTPUT_PATH ${${PROJECT_NAME}_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )
set(EXECUTABLE_OUTPUT_PATH ${${PROJECT_NAME}_BINARY_DIR}/bin CACHE PATH "Output directory for applications" )

# Group projects in "folders"
# ===================================================
if (CMAKE_VERSION VERSION_GREATER "2.8.0")
	set(ENABLE_SOLUTION_FOLDERS ON CACHE BOOL "Group projects under virtual folders of the compiler IDE (e.g. VisualStudio)")
endif (CMAKE_VERSION VERSION_GREATER "2.8.0")

if (ENABLE_SOLUTION_FOLDERS)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakeTargets")
endif (ENABLE_SOLUTION_FOLDERS)


# --------------------------------------------------------------------------
#   The list of "libs" which can be included can be found in:
#     http://www.mrpt.org/Libraries
#
#   The dependencies of a library are automatically added, so you only
#    need to specify the top-most libraries your code depend on.
# --------------------------------------------------------------------------
find_package( MRPT REQUIRED gui opengl base)

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif(CMAKE_COMPILER_IS_GNUCXX)

# Set optimized building:
if(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")


# ------------------------------------------------------------------
# Detect SuiteSparse libraries:
# If not found automatically, set SuiteSparse_DIR in CMake to the
# directory where SuiteSparse was built.
# ------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake") # Add the directory where FindSuiteSparse.cmake module can be found.

set(SuiteSparse_USE_LAPACK_BLAS ON)
find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
if(NOT SuiteSparse_FOUND)
        #set(SuiteSparse_VERBOSE ON)
        find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
else()
        message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
        include(${USE_SuiteSparse})
endif()
#message(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")
# ------------------------------------------------------------------
#   End of SuiteSparse detection
# ------------------------------------------------------------------

# We want libraries to be named "libXXX" and "libXXXd" in all compilers:
# ------------------------------------------------------------------------
set(CMAKE_DEBUG_POSTFIX  "d")
if(MSVC)
	set(SP_LIB_PREFIX "lib")  # Libs are: "libXXX"
endif(MSVC)

# LIB : libsparsembs
# --------------------------------
add_subdirectory(libsparsembs)

# Examples:
# --------------------------------
add_subdirectory(examples)
