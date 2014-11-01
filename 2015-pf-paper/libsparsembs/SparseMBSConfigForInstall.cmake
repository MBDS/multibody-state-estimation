##  SparseMBSConfigForInstall.cmake.in
##  Based on file by jesnault (jerome.esnault@inria.fr). 
##  Modified by Jose Luis Blanco (UAL)
##
##  File which define the USE_SparseMBS cmake variable for another project.
##  When cmake export a project cmake will generate an export file containing all the project's targets imported.
##  Here we just set the cmake variable (USE_SparseMBS) pointing to the file which do the import SparseMBS project stuff.
##
##  Usage for project which try to use SparseMBS :
##  (1) Set your SparseMBS_DIR to the dir containing this file.
##  (2) Then, in your CMakeLists.txt
##      find_package(SparseMBS NO_MODULE) ## NO_MODULE is optional (to bypass the FindSparseMBS.cmake is exist)
##      include(${USE_SparseMBS}) ## see UseSparseMBS.cmake for more infos (it does the include_directories)
##	(3) Then, in your target project you can use the SparseMBS_LIBRARIES
##

get_filename_component(SparseMBS_IMPORT_PREFIX 	"${CMAKE_CURRENT_LIST_FILE}" PATH)

## check which build system version we have to load (32 or 64 bits)
if(CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(SparseMBS_LIB_POSTFIX "64")## suffix for 32/64 inst dir placement
else()
  set(SparseMBS_LIB_POSTFIX "" ) ## suffix for 32/64 inst dir placement
endif()

set(USE_SparseMBS ${SparseMBS_IMPORT_PREFIX}/UseSparseMBS${SparseMBS_LIB_POSTFIX}.cmake)

if(EXISTS ${USE_SparseMBS})
	## do nothing, it's OK
else()
	message(SEND_ERROR "correct version of SparseMBS not found :\nUSE_SparseMBS=${USE_SparseMBS}")
	set(SparseMBS_FOUND OFF)
	set(SparseMBS_FOUND OFF)
endif()