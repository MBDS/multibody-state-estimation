PROJECT(libsparsembs)

## get POSTFIX for lib install dir
if(CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(LIB_POSTFIX "64" CACHE STRING "suffix for 32/64 inst dir placement")
else()
  set(LIB_POSTFIX "" CACHE STRING "suffix for 32/64 inst dir placement")
endif()
mark_as_advanced(LIB_POSTFIX)

include_directories("include")
SET(SPARSEMBS_INCLUDE_DIRS
	"${libsparsembs_SOURCE_DIR}/include"
	${SuiteSparse_INCLUDE_DIRS}
	CACHE INTERNAL "")

SET(SPARSEMBS_LIB_DIRS
	"${SuiteSparse_LIBS_DIR}"
	CACHE INTERNAL "")

# Tests
SET( BUILD_TESTING ON CACHE BOOL "Build unit tests")

# add before "D__MINGW32__"!!
IF( BUILD_TESTING)
	add_subdirectory(tests)
ENDIF( BUILD_TESTING)


# Used even for MSVC, so blas function names have the final "_"
#if (MSVC)
#	ADD_DEFINITIONS(-D__MINGW32__)
#endif(MSVC)


FILE(GLOB LIBHDRS "src/*.h")
FILE(GLOB LIBSRCS "src/*.cpp")
FILE(GLOB LIBHDRSPUB "include/sparsembs/*.h")
ADD_LIBRARY(sparsembs ${LIBHDRS} ${LIBSRCS} ${LIBHDRSPUB})

set_target_properties(sparsembs PROPERTIES FOLDER "Libraries")


SET(DEPENCENCIES_SPARSEMBS
	sparsembs
	CACHE INTERNAL "")

ADD_DEPENDENCIES(sparsembs ${DEPENCENCIES_SPARSEMBS})

# List of linking libs:
SET(SPARSEMBS_LIBS
	""
	CACHE INTERNAL "")

SET(SPARSEMBS_LIBS
	${SPARSEMBS_LIBS}
	${SuiteSparse_LIBRARIES}
	CACHE INTERNAL "")

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" GREATER 2.8.12)
	cmake_policy(SET CMP0023 NEW) ## just for respecting the new target_link_libraries(...) signature procedure
endif()

IF(BUILD_SHARED)
	TARGET_LINK_LIBRARIES(sparsembs 
		LINK_PRIVATE ${MRPT_LIBS} ${SPARSEMBS_LIBS}
		)
ELSE()
	TARGET_LINK_LIBRARIES(sparsembs 
		LINK_PUBLIC ${MRPT_LIBS} ${SPARSEMBS_LIBS}
		)
ENDIF()

SET(SPARSEMBS_LIBS
	sparsembs
	${SPARSEMBS_LIBS}
	CACHE INTERNAL "")

	
# Install target:
set_target_properties(sparsembs PROPERTIES PUBLIC_HEADER "${LIBHDRSPUB}")
install(TARGETS	sparsembs
		EXPORT 	sparsembs
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib${LIB_POSTFIX}
		ARCHIVE DESTINATION lib${LIB_POSTFIX}
		PUBLIC_HEADER DESTINATION include/sparsembs
)

## Install *Config files:
install(
	FILES 		SparseMBSConfigForInstall.cmake
	DESTINATION	.
	RENAME 		SparseMBSConfig.cmake 
	)
configure_file(UseSparseMBSForInstall.cmake.in UseSparseMBSForInstall.cmake @ONLY)
install(FILES 		"${libsparsembs_BINARY_DIR}/UseSparseMBSForInstall.cmake"
		DESTINATION	.
		RENAME 		UseSparseMBS${LIB_POSTFIX}.cmake
)

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" GREATER 2.8.12)
	cmake_policy(SET CMP0022 NEW) ## to not export the link of private targets
	set(additionalExportOption EXPORT_LINK_INTERFACE_LIBRARIES)
else()
	set(additionalExportOption )
endif()
## do the EXPORT for allowing other project to easily use this lib with cmake
install(EXPORT  	sparsembs
        DESTINATION cmake
        FILE        SparseMBS-config${LIB_POSTFIX}.cmake
		${additionalExportOption}
)
