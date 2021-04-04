project(mbse-fg-smoother-forward-dynamics)

find_package(mrpt-tclap)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mbse::mbse mrpt::tclap mrpt::system)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Apps")
install(TARGETS ${PROJECT_NAME})
