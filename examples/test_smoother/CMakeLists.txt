project(test_smoother)

include_directories(${SPARSEMBS_INCLUDE_DIRS})
link_directories(${SPARSEMBS_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mbse::mbse)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Examples")
