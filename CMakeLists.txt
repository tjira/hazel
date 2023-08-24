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

# add build target before all tests
add_test(build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})

# test: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b sto-3g -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf)
set_tests_properties(water_0-1_sto-3g_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.96590121")

# test: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2 COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b sto-3g -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf mp2)
set_tests_properties(water_0-1_sto-3g_mp2 PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.00485494")

# test: water 0 1 6-31G* HF
add_test(NAME water_0-1_6-31gs_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b 6-31g* -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf)
set_tests_properties(water_0-1_6-31gs_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.00679980")

# test: water 0 1 6-31G* MP2
add_test(NAME water_0-1_6-31gs_mp2 COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b 6-31g* -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf mp2)
set_tests_properties(water_0-1_6-31gs_mp2 PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.19820599406358")

# test: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b cc-pvdz -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf)
set_tests_properties(water_0-1_cc-pvdz_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02312305")

# test: water 0 1 CC-PVDZ MP2
add_test(NAME water_0-1_cc-pvdz_mp2 COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b cc-pvdz -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 1 hf mp2)
set_tests_properties(water_0-1_cc-pvdz_mp2 PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.22985471")

# test: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b sto-3g -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 2 hf)
set_tests_properties(water_1-2_sto-3g_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.66228146")

# test: water 1 2 6-31G* HF
add_test(NAME water_1-2_6-31gs_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b 6-31g* -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 2 hf)
set_tests_properties(water_1-2_6-31gs_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61187809")

# test: water 1 2 CC-PVDZ HF
add_test(NAME water_1-2_cc-pvdz_hf COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b cc-pvdz -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -s 2 hf)
set_tests_properties(water_1-2_cc-pvdz_hf PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63182597")
