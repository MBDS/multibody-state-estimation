project(mbse LANGUAGES CXX VERSION 0.1.0)

file(GLOB_RECURSE LIBHDRS "src/*.h")
file(GLOB_RECURSE LIBSRCS "src/*.cpp")
file(GLOB_RECURSE LIBHDRSPUB "include/mbse/*.h")
add_library(${PROJECT_NAME} SHARED
    ${LIBHDRS}
    ${LIBSRCS}
    ${LIBHDRSPUB}
)

add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})


set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Libraries")

target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${${PROJECT_NAME}_SOURCE_DIR}/src
)

if (GTSAM_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC gtsam gtsam_unstable)
  target_compile_definitions(${PROJECT_NAME} PUBLIC SPARSEMBS_HAVE_GTSAM)
endif()

if (SuiteSparse_FOUND)
    target_compile_definitions(${PROJECT_NAME} PUBLIC SPARSEMBS_HAVE_SUITESPARSE)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${SuiteSparse_LIBRARIES})
    message(STATUS "SuiteSparse_LIBRARIES: ${SuiteSparse_LIBRARIES}")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC ${MRPT_LIBRARIES})

# Shared options between GCC and CLANG:
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
	target_compile_options(${PROJECT_NAME} PRIVATE
        -Wall -pedantic
    )

	if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_compile_options(${PROJECT_NAME} PRIVATE
            "-O3"
        )
	endif()
endif()

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE
        /wd4275
        /wd4251
   )
endif()

# Tests ===========
set(BUILD_TESTING ON CACHE BOOL "Build unit tests")

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# ==== Install & export target ========
install(TARGETS	${PROJECT_NAME}
		EXPORT 	${PROJECT_NAME}-targets
		RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ${PROJECT_NAME}-targets
    FILE ${PROJECT_NAME}-targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
export(
	TARGETS ${PROJECT_NAME}
	FILE "${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}-targets.cmake"
	NAMESPACE mbse::
)


include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in
    ${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}-config-version.cmake
    VERSION ${YART_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}-config.cmake
        ${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
