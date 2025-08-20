cmake_minimum_required(VERSION 3.7)

project(DG_demo)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)

add_subdirectory(1D_scalar)
add_subdirectory(1D_Euler)
add_subdirectory(2D_Euler)