project(pf_test1)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mbse::mbse)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Apps")
