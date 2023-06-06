cmake_minimum_required(VERSION 3.15)
project(integral)

# set cmake variables
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(integral src/integral.cpp)
