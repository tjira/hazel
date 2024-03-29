cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 20)
include(FetchContent)
project(Hazel)

# set policies
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
if (${CMAKE_VERSION} VERSION_GREATER "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()

# set cmake variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -fprofile-arcs -ftest-coverage -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O3")

# static link on windows
if (WIN32)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -static")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -static")
else()
    if (NOT DEFINED NOAVX OR NOT ${NOAVX})
        string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx")
    endif()
    unset(NOAVX CACHE)
endif()

# set the variable that contains the compile flags
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

# declare libraries
FetchContent_Declare(argparse SYSTEM GIT_REPOSITORY https://github.com/p-ranav/argparse.git GIT_TAG ac4c578f6020d5164f66ccb26b4727dea657a12b)
FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 6f9ad7da6122fdb4197c0b43dfec09ec3525305e)
FetchContent_Declare(json GIT_REPOSITORY https://github.com/nlohmann/json.git GIT_TAG a259ecc51e1951e12f757ce17db958e9881e9c6c)
FetchContent_Declare(libint URL https://github.com/evaleev/libint/releases/download/v2.8.1/libint-2.8.1-mpqc4.tgz)

# fetch the libraries
FetchContent_MakeAvailable(argparse json)

# download dependencies if needed
if (${STANDALONE})
    FetchContent_MakeAvailable(eigen libint)
    if (WIN32)
        file(READ ${libint_SOURCE_DIR}/include/libint2/boys.h CONTENT)
        string(REPLACE "#define posix_memfree(p) ((free((p))))" "#define posix_memfree(p)" CONTENT "${CONTENT}")
        file(WRITE ${libint_SOURCE_DIR}/include/libint2/boys.h "${CONTENT}")
    endif()
else()
    find_package(Eigen3 REQUIRED)
endif()

# include directories (the second include is there because of boost/preprocessor.hpp)
include_directories(include ${CMAKE_BINARY_DIR}/include ${argparse_SOURCE_DIR}/include/argparse ${json_SOURCE_DIR}/include)

# find OpenMP if enabled
if (${GOMP})
    find_package(OpenMP REQUIRED)
endif()
unset(GOMP CACHE)

# add hazel executable
add_executable(hazel
    src/bagel.cpp
    src/ci.cpp
    src/determinant.cpp
    src/distributor.cpp
    src/eigen.cpp
    src/gradient.cpp
    src/hessian.cpp
    src/hf.cpp
    src/integral.cpp
    src/lambda.cpp
    src/main.cpp
    src/md.cpp
    src/mp.cpp
    src/optimizer.cpp
    src/orca.cpp
    src/qd.cpp
    src/parser.cpp
    src/printer.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
    src/utility.cpp
    src/wavetool.cpp
)
target_compile_definitions(hazel PRIVATE CXXFLAGS="${CXXFLAGS}" DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(hazel Eigen3::Eigen)

# link libint
if (${STANDALONE})
    target_link_libraries(hazel libint2)
else()
    target_link_libraries(hazel int2.a)
endif()

# link OpenMP on build type
if (OPENMP_FOUND AND CMAKE_BUILD_TYPE STREQUAL Release)
    target_link_libraries(hazel OpenMP::OpenMP_CXX)
endif()

# set default tools flags
set(HPLOT ON CACHE BOOL "")
set(HVIEW ON CACHE BOOL "")

# include tools
if (${HPLOT})
    include(tool/hplot/CMakeLists.txt)
endif()
unset(HPLOT CACHE)
if (${HVIEW})
    include(tool/hview/CMakeLists.txt)
endif()
unset(HVIEW CACHE)

# enable testing
enable_testing()

# add build target before all tests
add_test(build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})

# include the file with tests
include(CMakeTests.txt)
