cmake_minimum_required(VERSION 3.12)

project(tomfem)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -pedantic -Werror)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_library(libtomfem tomfem.cpp solve_linear_system.cpp)
target_include_directories(libtomfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(libtomfem Eigen3::Eigen)

add_executable(tomfem main.cpp)
target_link_libraries(tomfem libtomfem)
