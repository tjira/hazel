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
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -static -static-libgcc -static-libstdc++")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -static -static-libgcc -static-libstdc++")
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
add_executable(hazel src/ci.cpp src/determinant.cpp src/distributor.cpp src/eigen.cpp src/gradient.cpp src/hessian.cpp src/hf.cpp src/integral.cpp src/lambda.cpp src/main.cpp src/md.cpp src/mp.cpp src/optimizer.cpp src/qd.cpp src/printer.cpp src/system.cpp src/timer.cpp src/transform.cpp src/utility.cpp)
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
set_tests_properties(water_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.46449227")

# test-energy: water 0 1 MINI MP2
add_test(NAME water_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.49944113")

# test-gradient: water 0 1 MINI HF
add_test(NAME water_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.20e-01")

# test-energy: water 1 2 MINI HF
add_test(NAME water_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.05104959")

# test-energy: water -1 2 MINI HF
add_test(NAME water_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.12985091")

# test-energy: water 2 1 MINI HF
add_test(NAME water_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.97492639")

# test-energy: water 0 1 MIDI HF
add_test(NAME water_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.55206314")

# test-energy: water 0 1 MIDI MP2
add_test(NAME water_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.67849541")

# test-gradient: water 0 1 MIDI HF
add_test(NAME water_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.94e-02")

# test-energy: water 1 2 MIDI HF
add_test(NAME water_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.16473292")

# test-energy: water -1 2 MIDI HF
add_test(NAME water_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.35783998")

# test-energy: water 2 1 MIDI HF
add_test(NAME water_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.98230998")

# test-energy: water 0 1 STO-2G HF
add_test(NAME water_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.74384491")

# test-energy: water 0 1 STO-2G MP2
add_test(NAME water_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -72.78196295")

# test-gradient: water 0 1 STO-2G HF
add_test(NAME water_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.26e-02")

# test-energy: water 1 2 STO-2G HF
add_test(NAME water_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.47602203")

# test-energy: water -1 2 STO-2G HF
add_test(NAME water_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -72.14218024")

# test-energy: water 2 1 STO-2G HF
add_test(NAME water_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -71.48252894")

# test-energy: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.96590121")

# test-energy: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.00485494")

# test-gradient: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.59e-06")

# test-energy: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.66228146")

# test-energy: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.30165233")

# test-energy: water 2 1 STO-3G HF
add_test(NAME water_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -73.63296871")

# test-energy: water 0 1 STO-4G HF
add_test(NAME water_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.49925786")

# test-energy: water 0 1 STO-4G MP2
add_test(NAME water_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.53849667")

# test-gradient: water 0 1 STO-4G HF
add_test(NAME water_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.73e-03")

# test-energy: water 1 2 STO-4G HF
add_test(NAME water_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19143619")

# test-energy: water -1 2 STO-4G HF
add_test(NAME water_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.84401523")

# test-energy: water 2 1 STO-4G HF
add_test(NAME water_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.15802904")

# test-energy: water 0 1 STO-5G HF
add_test(NAME water_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63894007")

# test-energy: water 0 1 STO-5G MP2
add_test(NAME water_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.67822261")

# test-gradient: water 0 1 STO-5G HF
add_test(NAME water_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.20e-03")

# test-energy: water 1 2 STO-5G HF
add_test(NAME water_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.32991627")

# test-energy: water -1 2 STO-5G HF
add_test(NAME water_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.08480840")

# test-energy: water 2 1 STO-5G HF
add_test(NAME water_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.29530718")

# test-energy: water 0 1 STO-6G HF
add_test(NAME water_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.68117984")

# test-energy: water 0 1 STO-6G MP2
add_test(NAME water_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.72047253")

# test-gradient: water 0 1 STO-6G HF
add_test(NAME water_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.45e-03")

# test-energy: water 1 2 STO-6G HF
add_test(NAME water_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.37191401")

# test-energy: water -1 2 STO-6G HF
add_test(NAME water_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.12747742")

# test-energy: water 2 1 STO-6G HF
add_test(NAME water_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.33708116")

# test-energy: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.58368698")

# test-energy: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.70882320")

# test-gradient: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.89e-02")

# test-energy: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19927942")

# test-energy: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.34401430")

# test-energy: water 2 1 3-21G HF
add_test(NAME water_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.01546119")

# test-energy: water 0 1 4-31G HF
add_test(NAME water_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.90332760")

# test-energy: water 0 1 4-31G MP2
add_test(NAME water_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.03604979")

# test-gradient: water 0 1 4-31G HF
add_test(NAME water_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.96e-02")

# test-energy: water 1 2 4-31G HF
add_test(NAME water_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.50370193")

# test-energy: water -1 2 4-31G HF
add_test(NAME water_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.71458647")

# test-energy: water 2 1 4-31G HF
add_test(NAME water_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.32184569")

# test-energy: water 0 1 6-21G HF
add_test(NAME water_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88636231")

# test-energy: water 0 1 6-21G MP2
add_test(NAME water_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.01125292")

# test-gradient: water 0 1 6-21G HF
add_test(NAME water_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.61e-02")

# test-energy: water 1 2 6-21G HF
add_test(NAME water_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.50384554")

# test-energy: water -1 2 6-21G HF
add_test(NAME water_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64696191")

# test-energy: water 2 1 6-21G HF
add_test(NAME water_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.32055098")

# test-energy: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97974730")

# test-energy: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.11188710")

# test-gradient: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.24e-02")

# test-energy: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57902812")

# test-energy: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.79551492")

# test-energy: water 2 1 6-31G HF
add_test(NAME water_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.39811656")

# test-energy: water 0 1 6-31+G HF
add_test(NAME water_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.98680631")

# test-energy: water 0 1 6-31+G MP2
add_test(NAME water_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.12539844")

# test-gradient: water 0 1 6-31+G HF
add_test(NAME water_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.51e-02")

# test-energy: water 1 2 6-31+G HF
add_test(NAME water_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57953130")

# test-energy: water -1 2 6-31+G HF
add_test(NAME water_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.85563164")

# test-energy: water 2 1 6-31+G HF
add_test(NAME water_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.39956677")

# test-energy: water 0 1 6-31++G HF
add_test(NAME water_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.98709883")

# test-energy: water 0 1 6-31++G MP2
add_test(NAME water_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.12575848")

# test-gradient: water 0 1 6-31++G HF
add_test(NAME water_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.40e-02")

# test-energy: water 1 2 6-31++G HF
add_test(NAME water_1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.58001811")

# test-energy: water -1 2 6-31++G HF
add_test(NAME water_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.95034658")

# test-energy: water 2 1 6-31++G HF
add_test(NAME water_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.40029580")

# test-energy: water 0 1 6-31+G* HF
add_test(NAME water_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01240610")

# test-energy: water 0 1 6-31+G* MP2
add_test(NAME water_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.20742205")

# test-gradient: water 0 1 6-31+G* HF
add_test(NAME water_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.44e-02")

# test-energy: water 1 2 6-31+G* HF
add_test(NAME water_1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61087401")

# test-energy: water -1 2 6-31+G* HF
add_test(NAME water_-1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.87334902")

# test-energy: water 2 1 6-31+G* HF
add_test(NAME water_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.41183803")

# test-energy: water 0 1 6-31++G* HF
add_test(NAME water_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01270322")

# test-energy: water 0 1 6-31++G* MP2
add_test(NAME water_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.20784955")

# test-gradient: water 0 1 6-31++G* HF
add_test(NAME water_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.36e-02")

# test-energy: water 1 2 6-31++G* HF
add_test(NAME water_1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61140365")

# test-energy: water -1 2 6-31++G* HF
add_test(NAME water_-1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.94089778")

# test-energy: water 2 1 6-31++G* HF
add_test(NAME water_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.41261248")

# test-energy: water 0 1 6-31+G** HF
add_test(NAME water_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02579698")

# test-energy: water 0 1 6-31+G** MP2
add_test(NAME water_0-1_6-31pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.23094325")

# test-gradient: water 0 1 6-31+G** HF
add_test(NAME water_0-1_6-31pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.49e-02")

# test-energy: water 1 2 6-31+G** HF
add_test(NAME water_1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.62437667")

# test-energy: water -1 2 6-31+G** HF
add_test(NAME water_-1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88659765")

# test-energy: water 2 1 6-31+G** HF
add_test(NAME water_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42891666")

# test-energy: water 0 1 6-31++G** HF
add_test(NAME water_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02595479")

# test-energy: water 0 1 6-31++G** MP2
add_test(NAME water_0-1_6-31ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.23123645")

# test-gradient: water 0 1 6-31++G** HF
add_test(NAME water_0-1_6-31ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.44e-02")

# test-energy: water 1 2 6-31++G** HF
add_test(NAME water_1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.62468805")

# test-energy: water -1 2 6-31++G** HF
add_test(NAME water_-1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.95433890")

# test-energy: water 2 1 6-31++G** HF
add_test(NAME water_2-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42935135")

# test-energy: water 0 1 6-311G HF
add_test(NAME water_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.00465981")

# test-energy: water 0 1 6-311G MP2
add_test(NAME water_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.16402968")

# test-gradient: water 0 1 6-311G HF
add_test(NAME water_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.88e-02")

# test-energy: water 1 2 6-311G HF
add_test(NAME water_1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60344636")

# test-energy: water -1 2 6-311G HF
add_test(NAME water_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.87513663")

# test-energy: water 2 1 6-311G HF
add_test(NAME water_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42396046")

# test-energy: water 0 1 6-311+G HF
add_test(NAME water_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01004307")

# test-energy: water 0 1 6-311+G MP2
add_test(NAME water_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.17436628")

# test-gradient: water 0 1 6-311+G HF
add_test(NAME water_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.87e-02")

# test-energy: water 1 2 6-311+G HF
add_test(NAME water_1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60397282")

# test-energy: water -1 2 6-311+G HF
add_test(NAME water_-1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.88978960")

# test-energy: water 2 1 6-311+G HF
add_test(NAME water_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42541786")

# test-energy: water 0 1 6-311++G HF
add_test(NAME water_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.01029317")

# test-energy: water 0 1 6-311++G MP2
add_test(NAME water_0-1_6-311ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.17469569")

# test-gradient: water 0 1 6-311++G HF
add_test(NAME water_0-1_6-311ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.85e-02")

# test-energy: water 1 2 6-311++G HF
add_test(NAME water_1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.60453953")

# test-energy: water -1 2 6-311++G HF
add_test(NAME water_-1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97278401")

# test-energy: water 2 1 6-311++G HF
add_test(NAME water_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.42624950")

# test-energy: water 0 1 6-311+G* HF
add_test(NAME water_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.03182250")

# test-energy: water 0 1 6-311+G* MP2
add_test(NAME water_0-1_6-311pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.26083285")

# test-gradient: water 0 1 6-311+G* HF
add_test(NAME water_0-1_6-311pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.98e-02")

# test-energy: water 1 2 6-311+G* HF
add_test(NAME water_1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63186967")

# test-energy: water -1 2 6-311+G* HF
add_test(NAME water_-1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.83121424")

# test-energy: water 2 1 6-311+G* HF
add_test(NAME water_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43660218")

# test-energy: water 0 1 6-311++G* HF
add_test(NAME water_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.03201893")

# test-energy: water 0 1 6-311++G* MP2
add_test(NAME water_0-1_6-311ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.26112830")

# test-gradient: water 0 1 6-311++G* HF
add_test(NAME water_0-1_6-311ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.99e-02")

# test-energy: water 1 2 6-311++G* HF
add_test(NAME water_1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63234823")

# test-energy: water -1 2 6-311++G* HF
add_test(NAME water_-1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.96043569")

# test-energy: water 2 1 6-311++G* HF
add_test(NAME water_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43735855")

# test-energy: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04815693")

# test-energy: water 0 1 6-311+G** MP2
add_test(NAME water_0-1_6-311pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.29193795")

# test-gradient: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.60e-02")

# test-energy: water 1 2 6-311+G** HF
add_test(NAME water_1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64712237")

# test-energy: water -1 2 6-311+G** HF
add_test(NAME water_-1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.84554659")

# test-energy: water 2 1 6-311+G** HF
add_test(NAME water_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45108756")

# test-energy: water 0 1 6-311++G** HF
add_test(NAME water_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04827538")

# test-energy: water 0 1 6-311++G** MP2
add_test(NAME water_0-1_6-311ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-311ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.29216427")

# test-gradient: water 0 1 6-311++G** HF
add_test(NAME water_0-1_6-311ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.61e-02")

# test-energy: water 1 2 6-311++G** HF
add_test(NAME water_1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.64745286")

# test-energy: water -1 2 6-311++G** HF
add_test(NAME water_-1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97642811")

# test-energy: water 2 1 6-311++G** HF
add_test(NAME water_2-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45164137")

# test-energy: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02312305")

# test-energy: water 0 1 CC-PVDZ MP2
add_test(NAME water_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.22985471")

# test-gradient: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_cc-pvdz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.49e-02")

# test-energy: water 1 2 CC-PVDZ HF
add_test(NAME water_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63182597")

# test-energy: water -1 2 CC-PVDZ HF
add_test(NAME water_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.85263322")

# test-energy: water 2 1 CC-PVDZ HF
add_test(NAME water_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43265405")

# test-energy: water 0 1 CC-PVTZ HF
add_test(NAME water_0-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.05261153")

# test-energy: water 0 1 CC-PVTZ MP2
add_test(NAME water_0-1_cc-pvtz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_cc-pvtz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.33034247")

# test-gradient: water 0 1 CC-PVTZ HF
add_test(NAME water_0-1_cc-pvtz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_cc-pvtz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.58e-02")

# test-energy: water 1 2 CC-PVTZ HF
add_test(NAME water_1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.65543561")

# test-energy: water -1 2 CC-PVTZ HF
add_test(NAME water_-1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.92102029")

# test-energy: water 2 1 CC-PVTZ HF
add_test(NAME water_2-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45984337")
