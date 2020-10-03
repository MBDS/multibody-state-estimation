#  +-------------------------------------------------------------------------+
#  |            Multi Body State Estimation (mbse) C++ library               |
#  |                                                                         |
#  | Copyright (C) 2014-2020 University of Almeria                           |
#  | Copyright (C) 2020 University of Salento                                |
#  | See README for list of authors and papers                               |
#  | Distributed under 3-clause BSD license                                  |
#  |  See: <https://opensource.org/licenses/BSD-3-Clause>                    |
#  +-------------------------------------------------------------------------+
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW) # Allow project(xxx VERSION a.b.c)
endif()

project(MBSE)

cmake_minimum_required(VERSION 3.0)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

set(LIBRARY_OUTPUT_PATH ${${PROJECT_NAME}_BINARY_DIR}/lib CACHE PATH "Output directory for libraries" )
set(EXECUTABLE_OUTPUT_PATH ${${PROJECT_NAME}_BINARY_DIR}/bin CACHE PATH "Output directory for applications" )

# Group projects in "folders"
# ===================================================
set(ENABLE_SOLUTION_FOLDERS ON CACHE BOOL "Group projects under virtual folders of the compiler IDE (e.g. VisualStudio)")

if (ENABLE_SOLUTION_FOLDERS)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakeTargets")
endif ()

find_package(MRPT REQUIRED gui opengl)
find_package(GTSAM REQUIRED)


# Set optimized building:
if(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

# ------------------------------------------------------------------
# Detect SuiteSparse libraries:
# If not found automatically, set SuiteSparse_DIR in CMake to the
# directory where SuiteSparse was built.
# ------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake") # Add the directory where FindSuiteSparse.cmake module can be found.

option(DISABLE_SUITESPARSE OFF)
if (NOT DISABLE_SUITESPARSE)
    set(SuiteSparse_USE_LAPACK_BLAS ON)
    find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
    if(NOT SuiteSparse_FOUND)
        #set(SuiteSparse_VERBOSE ON)
        find_package(SuiteSparse QUIET) #REQUIRED # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
    endif()
    #message(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")
endif()

# Lib postfix for debug builds:
set(CMAKE_DEBUG_POSTFIX  "-dbg")
if(MSVC)
	set(CMAKE_SHARED_LIBRARY_PREFIX "lib")  # Libs are: "libXXX"
endif()
# To generate DLLs:
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

enable_testing()

# LIB : MBSE
# --------------------------------
add_subdirectory(libmbse)

# Examples:
# --------------------------------
option(BUILD_EXAMPLES ON)
if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()


# Show summary:
# --------------------------------
message(STATUS "${PROJECT_NAME} configuration:")
message(STATUS " MRPT version       : ${MRPT_VERSION}")
message(STATUS " GTSAM version      : ${GTSAM_VERSION}")
message(STATUS " SuiteSparse_FOUND  : ${SuiteSparse_FOUND}")
