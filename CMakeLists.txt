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
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

# find packages and include directories
include_directories(include lib)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# add hazel executable
add_executable(hazel src/ci.cpp src/distributor.cpp src/dynamics.cpp src/eigen.cpp src/gradient.cpp src/hessian.cpp src/hf.cpp src/integral.cpp src/main.cpp src/mp.cpp src/optimizer.cpp src/qdyn.cpp src/system.cpp src/timer.cpp src/transform.cpp)
target_compile_definitions(hazel PRIVATE CXXFLAGS="${CXXFLAGS}" DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(hazel Eigen3::Eigen int2.a xc.a)

# link OpenMP on build type
if (CMAKE_BUILD_TYPE STREQUAL Release AND OPENMP_FOUND)
    target_link_libraries(hazel OpenMP::OpenMP_CXX)
endif()

# enable testing
enable_testing()

# get all the test source files and create the source lists
file(GLOB SRC_TEST RELATIVE ${PROJECT_SOURCE_DIR} test/*.cpp)
create_test_sourcelist(TEST_LIST test_hazel.cpp ${SRC_TEST})

# add executables and set their properties
add_executable(test_hazel EXCLUDE_FROM_ALL src/dynamics.cpp src/eigen.cpp src/gradient.cpp src/hessian.cpp src/hf.cpp src/integral.cpp src/mp.cpp src/optimizer.cpp src/system.cpp src/transform.cpp src/timer.cpp ${TEST_LIST})
set_target_properties(test_hazel PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/test)
target_compile_definitions(test_hazel PRIVATE DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(test_hazel Eigen3::Eigen int2.a xc.a)

# extract energy test name
foreach (FILE ${SRC_TEST})
    get_filename_component(TNAME ${FILE} NAME_WE)
    list(APPEND TESTS test/${TNAME})
endforeach()

# add energy tests
foreach (TEST ${TESTS})
    add_test(NAME ${TEST} COMMAND test_hazel ${TEST})
endforeach()
