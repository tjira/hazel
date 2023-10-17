cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)
include(FetchContent)
project(Hazel)

# set policies
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
cmake_policy(SET CMP0135 NEW)

# set cmake variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_TESTING OFF)

# set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Og -fprofile-arcs -ftest-coverage -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O2")

# enable AVX
if (${AVX})
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -mavx")
endif()
unset(AVX CACHE)

# static link on windows
if (WIN32)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -static")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -static")
endif()

# set the variable that contains the compile flags
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

# download dependencies if needed
if (${STANDALONE})
    FetchContent_Declare(libint URL https://github.com/evaleev/libint/releases/download/v2.7.2/libint-2.7.2-mpqc4.tgz)
    FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4)
    FetchContent_MakeAvailable(eigen libint)
    if (WIN32)
        file(READ ${libint_SOURCE_DIR}/include/libint2/boys.h CONTENT)
        string(REPLACE "free(c);" "// free(c);" CONTENT "${CONTENT}")
        string(REPLACE "32ul)" "32ull)" CONTENT "${CONTENT}")
        file(WRITE ${libint_SOURCE_DIR}/include/libint2/boys.h "${CONTENT}")
    endif()
else()
    find_package(Eigen3 REQUIRED)
endif()

# include directories
include_directories(${CMAKE_BINARY_DIR}/include include lib)

# find OpenMP if enabled
if (${GOMP})
    find_package(OpenMP REQUIRED)
endif()
unset(GOMP CACHE)

# add hazel executable
add_executable(hazel
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
    src/qd.cpp
    src/printer.cpp
    src/system.cpp
    src/timer.cpp
    src/transform.cpp
    src/utility.cpp
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

# enable testing
enable_testing()

# add build target before all tests
add_test(build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})

# TESTS ==============================================================================================================================================

# test-energy: water 0 1 MINI HF
add_test(NAME water_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.46444890")

# test-energy: water 0 1 MINI MP2
add_test(NAME water_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.49939843")

# test-gradient: water 0 1 MINI HF
add_test(NAME water_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.20e-01")

# test-energy: water 1 2 MINI HF
add_test(NAME water_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.05094153")

# test-energy: water -1 2 MINI HF
add_test(NAME water_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.12970767")

# test-energy: water 2 1 MINI HF
add_test(NAME water_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.97473097")

# test-energy: water 0 1 MIDI HF
add_test(NAME water_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.55204757")

# test-energy: water 0 1 MIDI MP2
add_test(NAME water_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.67847718")

# test-gradient: water 0 1 MIDI HF
add_test(NAME water_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.94e-02")

# test-energy: water 1 2 MIDI HF
add_test(NAME water_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.16467325")

# test-energy: water -1 2 MIDI HF
add_test(NAME water_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.35780757")

# test-energy: water 2 1 MIDI HF
add_test(NAME water_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.98201205")

# test-energy: water 0 1 STO-2G HF
add_test(NAME water_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.74384824")

# test-energy: water 0 1 STO-2G MP2
add_test(NAME water_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -72.78196514")

# test-gradient: water 0 1 STO-2G HF
add_test(NAME water_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.27e-02")

# test-energy: water 1 2 STO-2G HF
add_test(NAME water_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.47596757")

# test-energy: water -1 2 STO-2G HF
add_test(NAME water_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.14198339")

# test-energy: water 2 1 STO-2G HF
add_test(NAME water_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -71.48238991")

# test-energy: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.96590106")

# test-energy: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.00485385")

# test-gradient: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.11e-04")

# test-energy: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.66222388")

# test-energy: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.30163897")

# test-energy: water 2 1 STO-3G HF
add_test(NAME water_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.63282683")

# test-energy: water 0 1 STO-4G HF
add_test(NAME water_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.49925928")

# test-energy: water 0 1 STO-4G MP2
add_test(NAME water_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.53849716")

# test-gradient: water 0 1 STO-4G HF
add_test(NAME water_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.56e-03")

# test-energy: water 1 2 STO-4G HF
add_test(NAME water_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19137932")

# test-energy: water -1 2 STO-4G HF
add_test(NAME water_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.84400568")

# test-energy: water 2 1 STO-4G HF
add_test(NAME water_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.15788713")

# test-energy: water 0 1 STO-5G HF
add_test(NAME water_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63894208")

# test-energy: water 0 1 STO-5G MP2
add_test(NAME water_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.67822366")

# test-gradient: water 0 1 STO-5G HF
add_test(NAME water_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.04e-03")

# test-energy: water 1 2 STO-5G HF
add_test(NAME water_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.32985988")

# test-energy: water -1 2 STO-5G HF
add_test(NAME water_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.08462494")

# test-energy: water 2 1 STO-5G HF
add_test(NAME water_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.29516555")

# test-energy: water 0 1 STO-6G HF
add_test(NAME water_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.68118164")

# test-energy: water 0 1 STO-6G MP2
add_test(NAME water_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.72047336")

# test-gradient: water 0 1 STO-6G HF
add_test(NAME water_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.29e-03")

# test-energy: water 1 2 STO-6G HF
add_test(NAME water_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.37185748")

# test-energy: water -1 2 STO-6G HF
add_test(NAME water_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.12729419")

# test-energy: water 2 1 STO-6G HF
add_test(NAME water_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.33693948")

# test-energy: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.58367190")

# test-energy: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.70880541")

# test-gradient: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.89e-02")

# test-energy: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19922162")

# test-energy: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.34397050")

# test-energy: water 2 1 3-21G HF
add_test(NAME water_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.01516788")

# test-energy: water 0 1 4-31G HF
add_test(NAME water_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.90331138")

# test-energy: water 0 1 4-31G MP2
add_test(NAME water_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.03602848")

# test-gradient: water 0 1 4-31G HF
add_test(NAME water_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.96e-02")

# test-energy: water 1 2 4-31G HF
add_test(NAME water_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.50364813")

# test-energy: water -1 2 4-31G HF
add_test(NAME water_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.71455526")

# test-energy: water 2 1 4-31G HF
add_test(NAME water_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.32155714")

# test-energy: water 0 1 6-21G HF
add_test(NAME water_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88634754")

# test-energy: water 0 1 6-21G MP2
add_test(NAME water_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.01123537")

# test-gradient: water 0 1 6-21G HF
add_test(NAME water_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.61e-02")

# test-energy: water 1 2 6-21G HF
add_test(NAME water_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.50378822")

# test-energy: water -1 2 6-21G HF
add_test(NAME water_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64691766")

# test-energy: water 2 1 6-21G HF
add_test(NAME water_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.32025587")

# test-energy: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97973091")

# test-energy: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.11186520")

# test-gradient: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.24e-02")

# test-energy: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57897498")

# test-energy: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.79548447")

# test-energy: water 2 1 6-31G HF
add_test(NAME water_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.39782841")

# test-energy: water 0 1 6-31+G HF
add_test(NAME water_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.98678744")

# test-energy: water 0 1 6-31+G MP2
add_test(NAME water_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.12537209")

# test-gradient: water 0 1 6-31+G HF
add_test(NAME water_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.51e-02")

# test-energy: water 1 2 6-31+G HF
add_test(NAME water_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57947879")

# test-energy: water -1 2 6-31+G HF
add_test(NAME water_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.85563394")

# test-energy: water 2 1 6-31+G HF
add_test(NAME water_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.39928044")

# test-energy: water 0 1 6-31++G HF
add_test(NAME water_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.98707984")

# test-energy: water 0 1 6-31++G MP2
add_test(NAME water_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.12573222")

# test-gradient: water 0 1 6-31++G HF
add_test(NAME water_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.39e-02")

# test-energy: water 1 2 6-31++G HF
add_test(NAME water_1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57996509")

# test-energy: water -1 2 6-31++G HF
add_test(NAME water_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.95033315")

# test-energy: water 2 1 6-31++G HF
add_test(NAME water_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.40000870")

# test-energy: water 0 1 6-31+G* HF
add_test(NAME water_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01240637")

# test-energy: water 0 1 6-31+G* MP2
add_test(NAME water_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.20741433")

# test-gradient: water 0 1 6-31+G* HF
add_test(NAME water_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.44e-02")

# test-energy: water 1 2 6-31+G* HF
add_test(NAME water_1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61084293")

# test-energy: water -1 2 6-31+G* HF
add_test(NAME water_-1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.87336981")

# test-energy: water 2 1 6-31+G* HF
add_test(NAME water_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.41155574")

# test-energy: water 0 1 6-31++G* HF
add_test(NAME water_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01270333")

# test-energy: water 0 1 6-31++G* MP2
add_test(NAME water_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.20784181")

# test-gradient: water 0 1 6-31++G* HF
add_test(NAME water_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.35e-02")

# test-energy: water 1 2 6-31++G* HF
add_test(NAME water_1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61137216")

# test-energy: water -1 2 6-31++G* HF
add_test(NAME water_-1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.94089222")

# test-energy: water 2 1 6-31++G* HF
add_test(NAME water_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.41232950")

# test-energy: water 0 1 6-31+G** HF
add_test(NAME water_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02580033")

# test-energy: water 0 1 6-31+G** MP2
add_test(NAME water_0-1_6-31pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.23094147")

# test-gradient: water 0 1 6-31+G** HF
add_test(NAME water_0-1_6-31pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.49e-02")

# test-energy: water 1 2 6-31+G** HF
add_test(NAME water_1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.62434789")

# test-energy: water -1 2 6-31+G** HF
add_test(NAME water_-1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88662199")

# test-energy: water 2 1 6-31+G** HF
add_test(NAME water_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42863636")

# test-energy: water 0 1 6-31++G** HF
add_test(NAME water_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02595812")

# test-energy: water 0 1 6-31++G** MP2
add_test(NAME water_0-1_6-31ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.23123478")

# test-gradient: water 0 1 6-31++G** HF
add_test(NAME water_0-1_6-31ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.44e-02")

# test-energy: water 1 2 6-31++G** HF
add_test(NAME water_1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.62465906")

# test-energy: water -1 2 6-31++G** HF
add_test(NAME water_-1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.95433626")

# test-energy: water 2 1 6-31++G** HF
add_test(NAME water_2-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42907071")

# test-energy: water 0 1 6-311G HF
add_test(NAME water_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.00464501")

# test-energy: water 0 1 6-311G MP2
add_test(NAME water_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.16400948")

# test-gradient: water 0 1 6-311G HF
add_test(NAME water_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.87e-02")

# test-energy: water 1 2 6-311G HF
add_test(NAME water_1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60339534")

# test-energy: water -1 2 6-311G HF
add_test(NAME water_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.87512375")

# test-energy: water 2 1 6-311G HF
add_test(NAME water_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42367641")

# test-energy: water 0 1 6-311+G HF
add_test(NAME water_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01002635")

# test-energy: water 0 1 6-311+G MP2
add_test(NAME water_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.17434251")

# test-gradient: water 0 1 6-311+G HF
add_test(NAME water_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.86e-02")

# test-energy: water 1 2 6-311+G HF
add_test(NAME water_1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60392195")

# test-energy: water -1 2 6-311+G HF
add_test(NAME water_-1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88977105")

# test-energy: water 2 1 6-311+G HF
add_test(NAME water_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42513386")

# test-energy: water 0 1 6-311++G HF
add_test(NAME water_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01027641")

# test-energy: water 0 1 6-311++G MP2
add_test(NAME water_0-1_6-311ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.17467187")

# test-gradient: water 0 1 6-311++G HF
add_test(NAME water_0-1_6-311ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.85e-02")

# test-energy: water 1 2 6-311++G HF
add_test(NAME water_1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60448858")

# test-energy: water -1 2 6-311++G HF
add_test(NAME water_-1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97277238")

# test-energy: water 2 1 6-311++G HF
add_test(NAME water_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42596524")

# test-energy: water 0 1 6-311+G* HF
add_test(NAME water_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.03182159")

# test-energy: water 0 1 6-311+G* MP2
add_test(NAME water_0-1_6-311pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.26082499")

# test-gradient: water 0 1 6-311+G* HF
add_test(NAME water_0-1_6-311pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.98e-02")

# test-energy: water 1 2 6-311+G* HF
add_test(NAME water_1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63183816")

# test-energy: water -1 2 6-311+G* HF
add_test(NAME water_-1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.83119117")

# test-energy: water 2 1 6-311+G* HF
add_test(NAME water_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43632169")

# test-energy: water 0 1 6-311++G* HF
add_test(NAME water_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.03201799")

# test-energy: water 0 1 6-311++G* MP2
add_test(NAME water_0-1_6-311ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.26112040")

# test-gradient: water 0 1 6-311++G* HF
add_test(NAME water_0-1_6-311ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.98e-02")

# test-energy: water 1 2 6-311++G* HF
add_test(NAME water_1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63231666")

# test-energy: water -1 2 6-311++G* HF
add_test(NAME water_-1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.96042982")

# test-energy: water 2 1 6-311++G* HF
add_test(NAME water_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43707787")

# test-energy: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04816368")

# test-energy: water 0 1 6-311+G** MP2
add_test(NAME water_0-1_6-311pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.29194340")

# test-gradient: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.59e-02")

# test-energy: water 1 2 6-311+G** HF
add_test(NAME water_1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64709681")

# test-energy: water -1 2 6-311+G** HF
add_test(NAME water_-1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.84552976")

# test-energy: water 2 1 6-311+G** HF
add_test(NAME water_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45080773")

# test-energy: water 0 1 6-311++G** HF
add_test(NAME water_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04828214")

# test-energy: water 0 1 6-311++G** MP2
add_test(NAME water_0-1_6-311ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.29216975")

# test-gradient: water 0 1 6-311++G** HF
add_test(NAME water_0-1_6-311ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.60e-02")

# test-energy: water 1 2 6-311++G** HF
add_test(NAME water_1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64742734")

# test-energy: water -1 2 6-311++G** HF
add_test(NAME water_-1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97642922")

# test-energy: water 2 1 6-311++G** HF
add_test(NAME water_2-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45136152")

# test-energy: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02313132")

# test-energy: water 0 1 CC-PVDZ MP2
add_test(NAME water_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.22986173")

# test-gradient: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_cc-pvdz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.48e-02")

# test-energy: water 1 2 CC-PVDZ HF
add_test(NAME water_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63180080")

# test-energy: water -1 2 CC-PVDZ HF
add_test(NAME water_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.85263131")

# test-energy: water 2 1 CC-PVDZ HF
add_test(NAME water_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43237168")

# test-energy: water 0 1 CC-PVTZ HF
add_test(NAME water_0-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.05261885")

# test-energy: water 0 1 CC-PVTZ MP2
add_test(NAME water_0-1_cc-pvtz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_cc-pvtz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.33034818")

# test-gradient: water 0 1 CC-PVTZ HF
add_test(NAME water_0-1_cc-pvtz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_cc-pvtz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.57e-02")

# test-energy: water 1 2 CC-PVTZ HF
add_test(NAME water_1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.65540917")

# test-energy: water -1 2 CC-PVTZ HF
add_test(NAME water_-1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.92102528")

# test-energy: water 2 1 CC-PVTZ HF
add_test(NAME water_2-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45956488")
