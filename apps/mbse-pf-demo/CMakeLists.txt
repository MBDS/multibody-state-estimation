project(mbse-pf-demo)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mbse::mbse)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Apps")
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
