cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)
project(Hazel)

# set cmake variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -fprofile-arcs -ftest-coverage -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O2")

# set the variable that contains the compile flags
if (${CMAKE_BUILD_TYPE} STREQUAL Debug)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
elseif (${CMAKE_BUILD_TYPE} STREQUAL Release)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

# find packages and include directories
include_directories(include lib)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# add hazel executable
add_executable(hazel src/distributor.cpp src/eigen.cpp src/integral.cpp src/main.cpp src/mp.cpp src/roothaan.cpp src/system.cpp src/timer.cpp src/transform.cpp)
target_compile_definitions(hazel PRIVATE CXXFLAGS="${CXXFLAGS}" DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(hazel Eigen3::Eigen int2.a)

# link OpenMP on build type
if (OPENMP_FOUND AND ${CMAKE_BUILD_TYPE} STREQUAL Release)
    target_link_libraries(hazel OpenMP::OpenMP_CXX)
endif()
