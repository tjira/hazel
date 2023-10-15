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

# test-energy: H2 0 1 MINI HF
add_test(NAME H2_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.07978768")

# test-energy: H2 0 1 MINI MP2
add_test(NAME H2_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.08897078")

# test-gradient: H2 0 1 MINI HF
add_test(NAME H2_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.74e-01")

# test-energy: H2 1 2 MINI HF
add_test(NAME H2_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.45406305")

# test-energy: H2 -1 2 MINI HF
add_test(NAME H2_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.66250229")

# test-energy: H2 2 1 MINI HF
add_test(NAME H2_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 MIDI HF
add_test(NAME H2_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12097505")

# test-energy: H2 0 1 MIDI MP2
add_test(NAME H2_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13754273")

# test-gradient: H2 0 1 MIDI HF
add_test(NAME H2_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.62e-02")

# test-energy: H2 1 2 MIDI HF
add_test(NAME H2_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54625147")

# test-energy: H2 -1 2 MIDI HF
add_test(NAME H2_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.89240862")

# test-energy: H2 2 1 MIDI HF
add_test(NAME H2_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 STO-2G HF
add_test(NAME H2_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.09375724")

# test-energy: H2 0 1 STO-2G MP2
add_test(NAME H2_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.10584673")

# test-gradient: H2 0 1 STO-2G HF
add_test(NAME H2_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.17e-02")

# test-energy: H2 1 2 STO-2G HF
add_test(NAME H2_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.51482974")

# test-energy: H2 -1 2 STO-2G HF
add_test(NAME H2_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.38158280")

# test-energy: H2 2 1 STO-2G HF
add_test(NAME H2_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 STO-3G HF
add_test(NAME H2_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.11750588")

# test-energy: H2 0 1 STO-3G MP2
add_test(NAME H2_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13000969")

# test-gradient: H2 0 1 STO-3G HF
add_test(NAME H2_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.00e-07")

# test-energy: H2 1 2 STO-3G HF
add_test(NAME H2_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.52728779")

# test-energy: H2 -1 2 STO-3G HF
add_test(NAME H2_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.41686011")

# test-energy: H2 2 1 STO-3G HF
add_test(NAME H2_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 STO-4G HF
add_test(NAME H2_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12396951")

# test-energy: H2 0 1 STO-4G MP2
add_test(NAME H2_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13648651")

# test-gradient: H2 0 1 STO-4G HF
add_test(NAME H2_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.83e-03")

# test-energy: H2 1 2 STO-4G HF
add_test(NAME H2_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.53048844")

# test-energy: H2 -1 2 STO-4G HF
add_test(NAME H2_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.42453067")

# test-energy: H2 2 1 STO-4G HF
add_test(NAME H2_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 STO-5G HF
add_test(NAME H2_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12565499")

# test-energy: H2 0 1 STO-5G MP2
add_test(NAME H2_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13818715")

# test-gradient: H2 0 1 STO-5G HF
add_test(NAME H2_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.34e-03")

# test-energy: H2 1 2 STO-5G HF
add_test(NAME H2_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.53134416")

# test-energy: H2 -1 2 STO-5G HF
add_test(NAME H2_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.42797462")

# test-energy: H2 2 1 STO-5G HF
add_test(NAME H2_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 STO-6G HF
add_test(NAME H2_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12621333")

# test-energy: H2 0 1 STO-6G MP2
add_test(NAME H2_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13874522")

# test-gradient: H2 0 1 STO-6G HF
add_test(NAME H2_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.62e-03")

# test-energy: H2 1 2 STO-6G HF
add_test(NAME H2_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.53162163")

# test-energy: H2 -1 2 STO-6G HF
add_test(NAME H2_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.42874104")

# test-energy: H2 2 1 STO-6G HF
add_test(NAME H2_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 3-21G HF
add_test(NAME H2_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12256431")

# test-energy: H2 0 1 3-21G MP2
add_test(NAME H2_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13961888")

# test-gradient: H2 0 1 3-21G HF
add_test(NAME H2_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.68e-02")

# test-energy: H2 1 2 3-21G HF
add_test(NAME H2_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54585789")

# test-energy: H2 -1 2 3-21G HF
add_test(NAME H2_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.85681864")

# test-energy: H2 2 1 3-21G HF
add_test(NAME H2_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 4-31G HF
add_test(NAME H2_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12658765")

# test-energy: H2 0 1 4-31G MP2
add_test(NAME H2_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14366746")

# test-gradient: H2 0 1 4-31G HF
add_test(NAME H2_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.07e-02")

# test-energy: H2 1 2 4-31G HF
add_test(NAME H2_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54888170")

# test-energy: H2 -1 2 4-31G HF
add_test(NAME H2_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.88714734")

# test-energy: H2 2 1 4-31G HF
add_test(NAME H2_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-21G HF
add_test(NAME H2_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12256431")

# test-energy: H2 0 1 6-21G MP2
add_test(NAME H2_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.13961888")

# test-gradient: H2 0 1 6-21G HF
add_test(NAME H2_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.68e-02")

# test-energy: H2 1 2 6-21G HF
add_test(NAME H2_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54585789")

# test-energy: H2 -1 2 6-21G HF
add_test(NAME H2_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.85681864")

# test-energy: H2 2 1 6-21G HF
add_test(NAME H2_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31G HF
add_test(NAME H2_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12658765")

# test-energy: H2 0 1 6-31G MP2
add_test(NAME H2_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14366746")

# test-gradient: H2 0 1 6-31G HF
add_test(NAME H2_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.07e-02")

# test-energy: H2 1 2 6-31G HF
add_test(NAME H2_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54888170")

# test-energy: H2 -1 2 6-31G HF
add_test(NAME H2_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.88714734")

# test-energy: H2 2 1 6-31G HF
add_test(NAME H2_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31+G HF
add_test(NAME H2_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12658765")

# test-energy: H2 0 1 6-31+G MP2
add_test(NAME H2_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14366746")

# test-gradient: H2 0 1 6-31+G HF
add_test(NAME H2_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.07e-02")

# test-energy: H2 1 2 6-31+G HF
add_test(NAME H2_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54888170")

# test-energy: H2 -1 2 6-31+G HF
add_test(NAME H2_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.88714734")

# test-energy: H2 2 1 6-31+G HF
add_test(NAME H2_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31++G HF
add_test(NAME H2_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12667628")

# test-energy: H2 0 1 6-31++G MP2
add_test(NAME H2_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14384299")

# test-gradient: H2 0 1 6-31++G HF
add_test(NAME H2_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.10e-02")

# test-energy: H2 1 2 6-31++G HF
add_test(NAME H2_1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54934470")

# test-energy: H2 -1 2 6-31++G HF
add_test(NAME H2_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.05095122")

# test-energy: H2 2 1 6-31++G HF
add_test(NAME H2_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31+G* HF
add_test(NAME H2_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12658765")

# test-energy: H2 0 1 6-31+G* MP2
add_test(NAME H2_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14366746")

# test-gradient: H2 0 1 6-31+G* HF
add_test(NAME H2_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.07e-02")

# test-energy: H2 1 2 6-31+G* HF
add_test(NAME H2_1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54888170")

# test-energy: H2 -1 2 6-31+G* HF
add_test(NAME H2_-1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.88714734")

# test-energy: H2 2 1 6-31+G* HF
add_test(NAME H2_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31++G* HF
add_test(NAME H2_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12667628")

# test-energy: H2 0 1 6-31++G* MP2
add_test(NAME H2_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14384299")

# test-gradient: H2 0 1 6-31++G* HF
add_test(NAME H2_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.10e-02")

# test-energy: H2 1 2 6-31++G* HF
add_test(NAME H2_1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.54934470")

# test-energy: H2 -1 2 6-31++G* HF
add_test(NAME H2_-1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.05095122")

# test-energy: H2 2 1 6-31++G* HF
add_test(NAME H2_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31+G** HF
add_test(NAME H2_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.13101744")

# test-energy: H2 0 1 6-31+G** MP2
add_test(NAME H2_0-1_6-31pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.15730711")

# test-gradient: H2 0 1 6-31+G** HF
add_test(NAME H2_0-1_6-31pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.38e-02")

# test-energy: H2 1 2 6-31+G** HF
add_test(NAME H2_1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55737959")

# test-energy: H2 -1 2 6-31+G** HF
add_test(NAME H2_-1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.89147489")

# test-energy: H2 2 1 6-31+G** HF
add_test(NAME H2_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-31++G** HF
add_test(NAME H2_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.13107518")

# test-energy: H2 0 1 6-31++G** MP2
add_test(NAME H2_0-1_6-31ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-31ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.15740571")

# test-gradient: H2 0 1 6-31++G** HF
add_test(NAME H2_0-1_6-31ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-31ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.42e-02")

# test-energy: H2 1 2 6-31++G** HF
add_test(NAME H2_1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55798629")

# test-energy: H2 -1 2 6-31++G** HF
add_test(NAME H2_-1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.05572985")

# test-energy: H2 2 1 6-31++G** HF
add_test(NAME H2_2-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311G HF
add_test(NAME H2_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12775899")

# test-energy: H2 0 1 6-311G MP2
add_test(NAME H2_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14543740")

# test-gradient: H2 0 1 6-311G HF
add_test(NAME H2_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.20e-02")

# test-energy: H2 1 2 6-311G HF
add_test(NAME H2_1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55043581")

# test-energy: H2 -1 2 6-311G HF
add_test(NAME H2_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.95928759")

# test-energy: H2 2 1 6-311G HF
add_test(NAME H2_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311+G HF
add_test(NAME H2_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12775899")

# test-energy: H2 0 1 6-311+G MP2
add_test(NAME H2_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14543740")

# test-gradient: H2 0 1 6-311+G HF
add_test(NAME H2_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.20e-02")

# test-energy: H2 1 2 6-311+G HF
add_test(NAME H2_1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55043581")

# test-energy: H2 -1 2 6-311+G HF
add_test(NAME H2_-1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.95928759")

# test-energy: H2 2 1 6-311+G HF
add_test(NAME H2_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311++G HF
add_test(NAME H2_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12775901")

# test-energy: H2 0 1 6-311++G MP2
add_test(NAME H2_0-1_6-311ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14546892")

# test-gradient: H2 0 1 6-311++G HF
add_test(NAME H2_0-1_6-311ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.20e-02")

# test-energy: H2 1 2 6-311++G HF
add_test(NAME H2_1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55050257")

# test-energy: H2 -1 2 6-311++G HF
add_test(NAME H2_-1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.05573854")

# test-energy: H2 2 1 6-311++G HF
add_test(NAME H2_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311+G* HF
add_test(NAME H2_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12775899")

# test-energy: H2 0 1 6-311+G* MP2
add_test(NAME H2_0-1_6-311pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14543740")

# test-gradient: H2 0 1 6-311+G* HF
add_test(NAME H2_0-1_6-311pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.20e-02")

# test-energy: H2 1 2 6-311+G* HF
add_test(NAME H2_1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55043581")

# test-energy: H2 -1 2 6-311+G* HF
add_test(NAME H2_-1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.95928759")

# test-energy: H2 2 1 6-311+G* HF
add_test(NAME H2_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311++G* HF
add_test(NAME H2_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12775901")

# test-energy: H2 0 1 6-311++G* MP2
add_test(NAME H2_0-1_6-311ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.14546892")

# test-gradient: H2 0 1 6-311++G* HF
add_test(NAME H2_0-1_6-311ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.20e-02")

# test-energy: H2 1 2 6-311++G* HF
add_test(NAME H2_1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55050257")

# test-energy: H2 -1 2 6-311++G* HF
add_test(NAME H2_-1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.05573854")

# test-energy: H2 2 1 6-311++G* HF
add_test(NAME H2_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311+G** HF
add_test(NAME H2_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.13208280")

# test-energy: H2 0 1 6-311+G** MP2
add_test(NAME H2_0-1_6-311pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.15976737")

# test-gradient: H2 0 1 6-311+G** HF
add_test(NAME H2_0-1_6-311pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.69e-02")

# test-energy: H2 1 2 6-311+G** HF
add_test(NAME H2_1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55983661")

# test-energy: H2 -1 2 6-311+G** HF
add_test(NAME H2_-1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.96667299")

# test-energy: H2 2 1 6-311+G** HF
add_test(NAME H2_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 6-311++G** HF
add_test(NAME H2_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.13210231")

# test-energy: H2 0 1 6-311++G** MP2
add_test(NAME H2_0-1_6-311ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_6-311ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.15979844")

# test-gradient: H2 0 1 6-311++G** HF
add_test(NAME H2_0-1_6-311ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_6-311ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.68e-02")

# test-energy: H2 1 2 6-311++G** HF
add_test(NAME H2_1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55984315")

# test-energy: H2 -1 2 6-311++G** HF
add_test(NAME H2_-1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.06037565")

# test-energy: H2 2 1 6-311++G** HF
add_test(NAME H2_2-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 CC-PVDZ HF
add_test(NAME H2_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.12776166")

# test-energy: H2 0 1 CC-PVDZ MP2
add_test(NAME H2_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.15387944")

# test-gradient: H2 0 1 CC-PVDZ HF
add_test(NAME H2_0-1_cc-pvdz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_cc-pvdz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.28e-02")

# test-energy: H2 1 2 CC-PVDZ HF
add_test(NAME H2_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.55601358")

# test-energy: H2 -1 2 CC-PVDZ HF
add_test(NAME H2_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.93129916")

# test-energy: H2 2 1 CC-PVDZ HF
add_test(NAME H2_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: H2 0 1 CC-PVTZ HF
add_test(NAME H2_0-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_0-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -1.13261973")

# test-energy: H2 0 1 CC-PVTZ MP2
add_test(NAME H2_0-1_cc-pvtz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(H2_0-1_cc-pvtz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -1.16419764")

# test-gradient: H2 0 1 CC-PVTZ HF
add_test(NAME H2_0-1_cc-pvtz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(H2_0-1_cc-pvtz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.55e-02")

# test-energy: H2 1 2 CC-PVTZ HF
add_test(NAME H2_1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.56071327")

# test-energy: H2 -1 2 CC-PVTZ HF
add_test(NAME H2_-1-2_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(H2_-1-2_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -0.96670697")

# test-energy: H2 2 1 CC-PVTZ HF
add_test(NAME H2_2-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/H2.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(H2_2-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: 0.74298642")

# test-energy: HF 0 1 MINI HF
add_test(NAME HF_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.35123900")

# test-energy: HF 0 1 MINI MP2
add_test(NAME HF_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.36557240")

# test-gradient: HF 0 1 MINI HF
add_test(NAME HF_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.09e-01")

# test-energy: HF 1 2 MINI HF
add_test(NAME HF_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.83559844")

# test-energy: HF -1 2 MINI HF
add_test(NAME HF_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.00335067")

# test-energy: HF 2 1 MINI HF
add_test(NAME HF_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.49386521")

# test-energy: HF 0 1 MIDI HF
add_test(NAME HF_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.42470921")

# test-energy: HF 0 1 MIDI MP2
add_test(NAME HF_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.55001918")

# test-gradient: HF 0 1 MIDI HF
add_test(NAME HF_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.00e-02")

# test-energy: HF 1 2 MIDI HF
add_test(NAME HF_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.92533490")

# test-energy: HF -1 2 MIDI HF
add_test(NAME HF_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.22428301")

# test-energy: HF 2 1 MIDI HF
add_test(NAME HF_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.43153060")

# test-energy: HF 0 1 STO-2G HF
add_test(NAME HF_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -95.60822105")

# test-energy: HF 0 1 STO-2G MP2
add_test(NAME HF_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -95.62746167")

# test-gradient: HF 0 1 STO-2G HF
add_test(NAME HF_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.86e-02")

# test-energy: HF 1 2 STO-2G HF
add_test(NAME HF_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -95.28534180")

# test-energy: HF -1 2 STO-2G HF
add_test(NAME HF_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -94.95864925")

# test-energy: HF 2 1 STO-2G HF
add_test(NAME HF_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -94.08998081")

# test-energy: HF 0 1 STO-3G HF
add_test(NAME HF_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.57284741")

# test-energy: HF 0 1 STO-3G MP2
add_test(NAME HF_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -98.59192832")

# test-gradient: HF 0 1 STO-3G HF
add_test(NAME HF_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.35e-06")

# test-energy: HF 1 2 STO-3G HF
add_test(NAME HF_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.18580665")

# test-energy: HF -1 2 STO-3G HF
add_test(NAME HF_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.98645033")

# test-energy: HF 2 1 STO-3G HF
add_test(NAME HF_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -96.92189540")

# test-energy: HF 0 1 STO-4G HF
add_test(NAME HF_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.26949853")

# test-energy: HF 0 1 STO-4G MP2
add_test(NAME HF_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.28860378")

# test-gradient: HF 0 1 STO-4G HF
add_test(NAME HF_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.70e-03")

# test-energy: HF 1 2 STO-4G HF
add_test(NAME HF_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.87476345")

# test-energy: HF -1 2 STO-4G HF
add_test(NAME HF_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.69233256")

# test-energy: HF 2 1 STO-4G HF
add_test(NAME HF_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.60254437")

# test-energy: HF 0 1 STO-5G HF
add_test(NAME HF_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.44808114")

# test-energy: HF 0 1 STO-5G MP2
add_test(NAME HF_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.46720153")

# test-gradient: HF 0 1 STO-5G HF
add_test(NAME HF_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.78e-03")

# test-energy: HF 1 2 STO-5G HF
add_test(NAME HF_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.05209640")

# test-energy: HF -1 2 STO-5G HF
add_test(NAME HF_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.87257465")

# test-energy: HF 2 1 STO-5G HF
add_test(NAME HF_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.77852391")

# test-energy: HF 0 1 STO-6G HF
add_test(NAME HF_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.50171677")

# test-energy: HF 0 1 STO-6G MP2
add_test(NAME HF_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.52084203")

# test-gradient: HF 0 1 STO-6G HF
add_test(NAME HF_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.65e-03")

# test-energy: HF 1 2 STO-6G HF
add_test(NAME HF_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.10544691")

# test-energy: HF -1 2 STO-6G HF
add_test(NAME HF_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.92660779")

# test-energy: HF 2 1 STO-6G HF
add_test(NAME HF_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.83161048")

# test-energy: HF 0 1 3-21G HF
add_test(NAME HF_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.45988643")

# test-energy: HF 0 1 3-21G MP2
add_test(NAME HF_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.58387089")

# test-gradient: HF 0 1 3-21G HF
add_test(NAME HF_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.71e-02")

# test-energy: HF 1 2 3-21G HF
add_test(NAME HF_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.96223019")

# test-energy: HF -1 2 3-21G HF
add_test(NAME HF_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.21694744")

# test-energy: HF 2 1 3-21G HF
add_test(NAME HF_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.46144391")

# test-energy: HF 0 1 4-31G HF
add_test(NAME HF_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.88616312")

# test-energy: HF 0 1 4-31G MP2
add_test(NAME HF_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.01756333")

# test-gradient: HF 0 1 4-31G HF
add_test(NAME HF_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.87e-02")

# test-energy: HF 1 2 4-31G HF
add_test(NAME HF_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.37032681")

# test-energy: HF -1 2 4-31G HF
add_test(NAME HF_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.69603810")

# test-energy: HF 2 1 4-31G HF
add_test(NAME HF_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.09961115")

# test-energy: HF 0 1 6-21G HF
add_test(NAME HF_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.84214315")

# test-energy: HF 0 1 6-21G MP2
add_test(NAME HF_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -99.96577749")

# test-gradient: HF 0 1 6-21G HF
add_test(NAME HF_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.38e-02")

# test-energy: HF 1 2 6-21G HF
add_test(NAME HF_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.34734714")

# test-energy: HF -1 2 6-21G HF
add_test(NAME HF_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.59862324")

# test-energy: HF 2 1 6-21G HF
add_test(NAME HF_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.84828719")

# test-energy: HF 0 1 6-31G HF
add_test(NAME HF_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.98220328")

# test-energy: HF 0 1 6-31G MP2
add_test(NAME HF_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.11293499")

# test-gradient: HF 0 1 6-31G HF
add_test(NAME HF_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.08e-02")

# test-energy: HF 1 2 6-31G HF
add_test(NAME HF_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.46447466")

# test-energy: HF -1 2 6-31G HF
add_test(NAME HF_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.79654708")

# test-energy: HF 2 1 6-31G HF
add_test(NAME HF_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.19505815")

# test-energy: HF 0 1 6-31+G HF
add_test(NAME HF_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.99442560")

# test-energy: HF 0 1 6-31+G MP2
add_test(NAME HF_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.13403331")

# test-gradient: HF 0 1 6-31+G HF
add_test(NAME HF_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.68e-02")

# test-energy: HF 1 2 6-31+G HF
add_test(NAME HF_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.46484943")

# test-energy: HF -1 2 6-31+G HF
add_test(NAME HF_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.82929218")

# test-energy: HF 2 1 6-31+G HF
add_test(NAME HF_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.96946528")

# test-energy: HF 0 1 6-31++G HF
add_test(NAME HF_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.99477530")

# test-energy: HF 0 1 6-31++G MP2
add_test(NAME HF_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.13436476")

# test-gradient: HF 0 1 6-31++G HF
add_test(NAME HF_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.68e-02")

# test-energy: HF 1 2 6-31++G HF
add_test(NAME HF_1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.46547475")

# test-energy: HF -1 2 6-31++G HF
add_test(NAME HF_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.95806923")

# test-energy: HF 2 1 6-31++G HF
add_test(NAME HF_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.97037186")

# test-energy: HF 0 1 6-31+G* HF
add_test(NAME HF_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01077673")

# test-energy: HF 0 1 6-31+G* MP2
add_test(NAME HF_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.20112775")

# test-gradient: HF 0 1 6-31+G* HF
add_test(NAME HF_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.28e-02")

# test-energy: HF 1 2 6-31+G* HF
add_test(NAME HF_1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.48528482")

# test-energy: HF -1 2 6-31+G* HF
add_test(NAME HF_-1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.83695680")

# test-energy: HF 2 1 6-31+G* HF
add_test(NAME HF_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.97736501")

# test-energy: HF 0 1 6-31++G* HF
add_test(NAME HF_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01106252")

# test-energy: HF 0 1 6-31++G* MP2
add_test(NAME HF_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.20142613")

# test-gradient: HF 0 1 6-31++G* HF
add_test(NAME HF_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.29e-02")

# test-energy: HF 1 2 6-31++G* HF
add_test(NAME HF_1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.48583408")

# test-energy: HF -1 2 6-31++G* HF
add_test(NAME HF_-1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.97163500")

# test-energy: HF 2 1 6-31++G* HF
add_test(NAME HF_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -97.97821574")

# test-energy: HF 0 1 6-31+G** HF
add_test(NAME HF_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01957264")

# test-energy: HF 0 1 6-31+G** MP2
add_test(NAME HF_0-1_6-31pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.21395831")

# test-gradient: HF 0 1 6-31+G** HF
add_test(NAME HF_0-1_6-31pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.60e-02")

# test-energy: HF 1 2 6-31+G** HF
add_test(NAME HF_1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.49462174")

# test-energy: HF -1 2 6-31+G** HF
add_test(NAME HF_-1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.84476182")

# test-energy: HF 2 1 6-31+G** HF
add_test(NAME HF_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.22606111")

# test-energy: HF 0 1 6-31++G** HF
add_test(NAME HF_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01970190")

# test-energy: HF 0 1 6-31++G** MP2
add_test(NAME HF_0-1_6-31ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-31ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.21410038")

# test-gradient: HF 0 1 6-31++G** HF
add_test(NAME HF_0-1_6-31ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-31ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.58e-02")

# test-energy: HF 1 2 6-31++G** HF
add_test(NAME HF_1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.49493745")

# test-energy: HF -1 2 6-31++G** HF
add_test(NAME HF_-1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.98044454")

# test-energy: HF 2 1 6-31++G** HF
add_test(NAME HF_2-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.22663317")

# test-energy: HF 0 1 6-311G HF
add_test(NAME HF_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01953854")

# test-energy: HF 0 1 6-311G MP2
add_test(NAME HF_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.18101911")

# test-gradient: HF 0 1 6-311G HF
add_test(NAME HF_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.59e-02")

# test-energy: HF 1 2 6-311G HF
add_test(NAME HF_1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.49729682")

# test-energy: HF -1 2 6-311G HF
add_test(NAME HF_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.89336428")

# test-energy: HF 2 1 6-311G HF
add_test(NAME HF_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.23023663")

# test-energy: HF 0 1 6-311+G HF
add_test(NAME HF_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.02606986")

# test-energy: HF 0 1 6-311+G MP2
add_test(NAME HF_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.19319424")

# test-gradient: HF 0 1 6-311+G HF
add_test(NAME HF_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.33e-02")

# test-energy: HF 1 2 6-311+G HF
add_test(NAME HF_1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.49903779")

# test-energy: HF -1 2 6-311+G HF
add_test(NAME HF_-1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.93103491")

# test-energy: HF 2 1 6-311+G HF
add_test(NAME HF_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.23453531")

# test-energy: HF 0 1 6-311++G HF
add_test(NAME HF_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.02610034")

# test-energy: HF 0 1 6-311++G MP2
add_test(NAME HF_0-1_6-311ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.19325728")

# test-gradient: HF 0 1 6-311++G HF
add_test(NAME HF_0-1_6-311ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.32e-02")

# test-energy: HF 1 2 6-311++G HF
add_test(NAME HF_1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.49910318")

# test-energy: HF -1 2 6-311++G HF
add_test(NAME HF_-1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.98904125")

# test-energy: HF 2 1 6-311++G HF
add_test(NAME HF_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.23462635")

# test-energy: HF 0 1 6-311+G* HF
add_test(NAME HF_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.03734148")

# test-energy: HF 0 1 6-311+G* MP2
add_test(NAME HF_0-1_6-311pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.27709881")

# test-gradient: HF 0 1 6-311+G* HF
add_test(NAME HF_0-1_6-311pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.59e-02")

# test-energy: HF 1 2 6-311+G* HF
add_test(NAME HF_1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.51487081")

# test-energy: HF -1 2 6-311+G* HF
add_test(NAME HF_-1-2_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.93855038")

# test-energy: HF 2 1 6-311+G* HF
add_test(NAME HF_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.24851137")

# test-energy: HF 0 1 6-311++G* HF
add_test(NAME HF_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.03736764")

# test-energy: HF 0 1 6-311++G* MP2
add_test(NAME HF_0-1_6-311ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.27716001")

# test-gradient: HF 0 1 6-311++G* HF
add_test(NAME HF_0-1_6-311ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.59e-02")

# test-energy: HF 1 2 6-311++G* HF
add_test(NAME HF_1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.51492162")

# test-energy: HF -1 2 6-311++G* HF
add_test(NAME HF_-1-2_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.99878897")

# test-energy: HF 2 1 6-311++G* HF
add_test(NAME HF_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.24858544")

# test-energy: HF 0 1 6-311+G** HF
add_test(NAME HF_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.04942503")

# test-energy: HF 0 1 6-311+G** MP2
add_test(NAME HF_0-1_6-311pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.29624105")

# test-gradient: HF 0 1 6-311+G** HF
add_test(NAME HF_0-1_6-311pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.22e-02")

# test-energy: HF 1 2 6-311+G** HF
add_test(NAME HF_1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.52625161")

# test-energy: HF -1 2 6-311+G** HF
add_test(NAME HF_-1-2_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.94893464")

# test-energy: HF 2 1 6-311+G** HF
add_test(NAME HF_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.25909550")

# test-energy: HF 0 1 6-311++G** HF
add_test(NAME HF_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.04947242")

# test-energy: HF 0 1 6-311++G** MP2
add_test(NAME HF_0-1_6-311ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_6-311ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.29634481")

# test-gradient: HF 0 1 6-311++G** HF
add_test(NAME HF_0-1_6-311ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_6-311ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.22e-02")

# test-energy: HF 1 2 6-311++G** HF
add_test(NAME HF_1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.52635316")

# test-energy: HF -1 2 6-311++G** HF
add_test(NAME HF_-1-2_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.00971767")

# test-energy: HF 2 1 6-311++G** HF
add_test(NAME HF_2-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.25923756")

# test-energy: HF 0 1 CC-PVDZ HF
add_test(NAME HF_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.01643971")

# test-energy: HF 0 1 CC-PVDZ MP2
add_test(NAME HF_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HF_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -100.22185388")

# test-gradient: HF 0 1 CC-PVDZ HF
add_test(NAME HF_0-1_cc-pvdz_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HF_0-1_cc-pvdz_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.49e-02")

# test-energy: HF 1 2 CC-PVDZ HF
add_test(NAME HF_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.50233471")

# test-energy: HF -1 2 CC-PVDZ HF
add_test(NAME HF_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HF_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -99.85059666")

# test-energy: HF 2 1 CC-PVDZ HF
add_test(NAME HF_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.23254202")

# test-energy: HF 0 1 CC-PVTZ HF
add_test(NAME HF_0-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_0-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -100.05470646")

# test-energy: HF 2 1 CC-PVTZ HF
add_test(NAME HF_2-1_cc-pvtz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvtz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HF.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HF_2-1_cc-pvtz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -98.26844468")

# test-energy: HCl 0 1 MINI HF
add_test(NAME HCl_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.83424600")

# test-energy: HCl 0 1 MINI MP2
add_test(NAME HCl_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -457.84534217")

# test-gradient: HCl 0 1 MINI HF
add_test(NAME HCl_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.11e-01")

# test-energy: HCl 1 2 MINI HF
add_test(NAME HCl_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.35826380")

# test-energy: HCl -1 2 MINI HF
add_test(NAME HCl_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.62030958")

# test-energy: HCl 2 1 MINI HF
add_test(NAME HCl_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -456.38106810")

# test-energy: HCl 0 1 MIDI HF
add_test(NAME HCl_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.88522797")

# test-energy: HCl 0 1 MIDI MP2
add_test(NAME HCl_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -457.90748625")

# test-gradient: HCl 0 1 MIDI HF
add_test(NAME HCl_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.82e-03")

# test-energy: HCl 1 2 MIDI HF
add_test(NAME HCl_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.43983761")

# test-energy: HCl -1 2 MIDI HF
add_test(NAME HCl_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.71475875")

# test-energy: HCl 2 1 MIDI HF
add_test(NAME HCl_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -456.50476131")

# test-energy: HCl 0 1 STO-2G HF
add_test(NAME HCl_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -442.79948377")

# test-energy: HCl 0 1 STO-2G MP2
add_test(NAME HCl_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -442.81269689")

# test-gradient: HCl 0 1 STO-2G HF
add_test(NAME HCl_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.16e-03")

# test-energy: HCl 1 2 STO-2G HF
add_test(NAME HCl_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -442.41531807")

# test-energy: HCl -1 2 STO-2G HF
add_test(NAME HCl_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -442.40819494")

# test-energy: HCl 2 1 STO-2G HF
add_test(NAME HCl_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -441.48469598")

# test-energy: HCl 0 1 STO-3G HF
add_test(NAME HCl_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -455.13601165")

# test-energy: HCl 0 1 STO-3G MP2
add_test(NAME HCl_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -455.14967272")

# test-gradient: HCl 0 1 STO-3G HF
add_test(NAME HCl_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.62e-06")

# test-energy: HCl 1 2 STO-3G HF
add_test(NAME HCl_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -454.73176041")

# test-energy: HCl -1 2 STO-3G HF
add_test(NAME HCl_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -454.76927930")

# test-energy: HCl 2 1 STO-3G HF
add_test(NAME HCl_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -453.77958828")

# test-energy: HCl 0 1 STO-4G HF
add_test(NAME HCl_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.70809459")

# test-energy: HCl 0 1 STO-4G MP2
add_test(NAME HCl_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -457.72185327")

# test-gradient: HCl 0 1 STO-4G HF
add_test(NAME HCl_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.07e-04")

# test-energy: HCl 1 2 STO-4G HF
add_test(NAME HCl_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.30182917")

# test-energy: HCl -1 2 STO-4G HF
add_test(NAME HCl_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.34518931")

# test-energy: HCl 2 1 STO-4G HF
add_test(NAME HCl_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -456.34706383")

# test-energy: HCl 0 1 STO-5G HF
add_test(NAME HCl_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.34960855")

# test-energy: HCl 0 1 STO-5G MP2
add_test(NAME HCl_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -458.36336260")

# test-gradient: HCl 0 1 STO-5G HF
add_test(NAME HCl_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.24e-03")

# test-energy: HCl 1 2 STO-5G HF
add_test(NAME HCl_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.94244604")

# test-energy: HCl -1 2 STO-5G HF
add_test(NAME HCl_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.98751928")

# test-energy: HCl 2 1 STO-5G HF
add_test(NAME HCl_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -456.98664425")

# test-energy: HCl 0 1 STO-6G HF
add_test(NAME HCl_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.54663321")

# test-energy: HCl 0 1 STO-6G MP2
add_test(NAME HCl_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -458.56037406")

# test-gradient: HCl 0 1 STO-6G HF
add_test(NAME HCl_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.38e-03")

# test-energy: HCl 1 2 STO-6G HF
add_test(NAME HCl_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.13862901")

# test-energy: HCl -1 2 STO-6G HF
add_test(NAME HCl_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.18543747")

# test-energy: HCl 2 1 STO-6G HF
add_test(NAME HCl_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.18194350")

# test-energy: HCl 0 1 3-21G HF
add_test(NAME HCl_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.86922732")

# test-energy: HCl 0 1 3-21G MP2
add_test(NAME HCl_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -457.92633131")

# test-gradient: HCl 0 1 3-21G HF
add_test(NAME HCl_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.51e-02")

# test-energy: HCl 1 2 3-21G HF
add_test(NAME HCl_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.43411618")

# test-energy: HCl -1 2 3-21G HF
add_test(NAME HCl_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -457.73077477")

# test-energy: HCl 2 1 3-21G HF
add_test(NAME HCl_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -456.52878326")

# test-energy: HCl 0 1 4-31G HF
add_test(NAME HCl_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.56350044")

# test-energy: HCl 0 1 4-31G MP2
add_test(NAME HCl_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -459.62218129")

# test-gradient: HCl 0 1 4-31G HF
add_test(NAME HCl_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.15e-02")

# test-energy: HCl 1 2 4-31G HF
add_test(NAME HCl_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.13298470")

# test-energy: HCl -1 2 4-31G HF
add_test(NAME HCl_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.42737770")

# test-energy: HCl 2 1 4-31G HF
add_test(NAME HCl_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.22433455")

# test-energy: HCl 0 1 6-21G HF
add_test(NAME HCl_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.02510840")

# test-energy: HCl 0 1 6-21G MP2
add_test(NAME HCl_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.08017971")

# test-gradient: HCl 0 1 6-21G HF
add_test(NAME HCl_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.50e-02")

# test-energy: HCl 1 2 6-21G HF
add_test(NAME HCl_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.59136398")

# test-energy: HCl -1 2 6-21G HF
add_test(NAME HCl_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.88513940")

# test-energy: HCl 2 1 6-21G HF
add_test(NAME HCl_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.68440699")

# test-energy: HCl 0 1 6-31G HF
add_test(NAME HCl_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.03699500")

# test-energy: HCl 0 1 6-31G MP2
add_test(NAME HCl_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.08845879")

# test-gradient: HCl 0 1 6-31G HF
add_test(NAME HCl_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.45e-02")

# test-energy: HCl 1 2 6-31G HF
add_test(NAME HCl_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.60460903")

# test-energy: HCl -1 2 6-31G HF
add_test(NAME HCl_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.91008676")

# test-energy: HCl 2 1 6-31G HF
add_test(NAME HCl_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.70053828")

# test-energy: HCl 0 1 6-31+G HF
add_test(NAME HCl_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.03824856")

# test-energy: HCl 0 1 6-31+G MP2
add_test(NAME HCl_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.09119101")

# test-gradient: HCl 0 1 6-31+G HF
add_test(NAME HCl_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.33e-02")

# test-energy: HCl 1 2 6-31+G HF
add_test(NAME HCl_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.60604785")

# test-energy: HCl -1 2 6-31+G HF
add_test(NAME HCl_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.96513012")

# test-energy: HCl 2 1 6-31+G HF
add_test(NAME HCl_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.70213465")

# test-energy: HCl 0 1 6-31++G HF
add_test(NAME HCl_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.03843546")

# test-energy: HCl 0 1 6-31++G MP2
add_test(NAME HCl_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.09151379")

# test-gradient: HCl 0 1 6-31++G HF
add_test(NAME HCl_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.29e-02")

# test-energy: HCl 1 2 6-31++G HF
add_test(NAME HCl_1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.60625059")

# test-energy: HCl -1 2 6-31++G HF
add_test(NAME HCl_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.00654261")

# test-energy: HCl 2 1 6-31++G HF
add_test(NAME HCl_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.70230627")

# test-energy: HCl 0 1 6-31+G* HF
add_test(NAME HCl_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.05815153")

# test-energy: HCl 0 1 6-31+G* MP2
add_test(NAME HCl_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.19703001")

# test-gradient: HCl 0 1 6-31+G* HF
add_test(NAME HCl_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.04e-02")

# test-energy: HCl 1 2 6-31+G* HF
add_test(NAME HCl_1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63339117")

# test-energy: HCl -1 2 6-31+G* HF
add_test(NAME HCl_-1-2_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.97899797")

# test-energy: HCl 2 1 6-31+G* HF
add_test(NAME HCl_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.72915874")

# test-energy: HCl 0 1 6-31++G* HF
add_test(NAME HCl_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.05835868")

# test-energy: HCl 0 1 6-31++G* MP2
add_test(NAME HCl_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.19737265")

# test-gradient: HCl 0 1 6-31++G* HF
add_test(NAME HCl_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.00e-02")

# test-energy: HCl 1 2 6-31++G* HF
add_test(NAME HCl_1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63362346")

# test-energy: HCl -1 2 6-31++G* HF
add_test(NAME HCl_-1-2_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.02311500")

# test-energy: HCl 2 1 6-31++G* HF
add_test(NAME HCl_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.72936534")

# test-energy: HCl 0 1 6-31+G** HF
add_test(NAME HCl_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.06458124")

# test-energy: HCl 0 1 6-31+G** MP2
add_test(NAME HCl_0-1_6-31pgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31pgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.21017119")

# test-gradient: HCl 0 1 6-31+G** HF
add_test(NAME HCl_0-1_6-31pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.19e-02")

# test-energy: HCl 1 2 6-31+G** HF
add_test(NAME HCl_1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63995635")

# test-energy: HCl -1 2 6-31+G** HF
add_test(NAME HCl_-1-2_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.98484833")

# test-energy: HCl 2 1 6-31+G** HF
add_test(NAME HCl_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.73591172")

# test-energy: HCl 0 1 6-31++G** HF
add_test(NAME HCl_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.06471081")

# test-energy: HCl 0 1 6-31++G** MP2
add_test(NAME HCl_0-1_6-31ppgss_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-31ppgss_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.21044654")

# test-gradient: HCl 0 1 6-31++G** HF
add_test(NAME HCl_0-1_6-31ppgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-31ppgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.18e-02")

# test-energy: HCl 1 2 6-31++G** HF
add_test(NAME HCl_1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.64010035")

# test-energy: HCl -1 2 6-31++G** HF
add_test(NAME HCl_-1-2_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.02947560")

# test-energy: HCl 2 1 6-31++G** HF
add_test(NAME HCl_2-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.73603339")

# test-energy: HCl 0 1 6-311G HF
add_test(NAME HCl_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.06657373")

# test-energy: HCl 0 1 6-311G MP2
add_test(NAME HCl_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.16910309")

# test-gradient: HCl 0 1 6-311G HF
add_test(NAME HCl_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.28e-03")

# test-energy: HCl 1 2 6-311G HF
add_test(NAME HCl_1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63379372")

# test-energy: HCl -1 2 6-311G HF
add_test(NAME HCl_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.97451543")

# test-energy: HCl 2 1 6-311G HF
add_test(NAME HCl_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.72944753")

# test-energy: HCl 0 1 6-311+G HF
add_test(NAME HCl_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.06745433")

# test-energy: HCl 0 1 6-311+G MP2
add_test(NAME HCl_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.17049448")

# test-gradient: HCl 0 1 6-311+G HF
add_test(NAME HCl_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.77e-03")

# test-energy: HCl 1 2 6-311+G HF
add_test(NAME HCl_1-2_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63566558")

# test-energy: HCl 2 1 6-311+G HF
add_test(NAME HCl_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.73103913")

# test-energy: HCl 0 1 6-311++G HF
add_test(NAME HCl_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.06770794")

# test-energy: HCl 0 1 6-311++G MP2
add_test(NAME HCl_0-1_6-311ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_6-311ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.17079636")

# test-gradient: HCl 0 1 6-311++G HF
add_test(NAME HCl_0-1_6-311ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-311ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.27e-03")

# test-energy: HCl 1 2 6-311++G HF
add_test(NAME HCl_1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.63606456")

# test-energy: HCl -1 2 6-311++G HF
add_test(NAME HCl_-1-2_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.03487390")

# test-energy: HCl 2 1 6-311++G HF
add_test(NAME HCl_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.73159626")

# test-energy: HCl 0 1 6-311+G* HF
add_test(NAME HCl_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.08612935")

# test-gradient: HCl 0 1 6-311+G* HF
add_test(NAME HCl_0-1_6-311pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-311pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.58e-02")

# test-energy: HCl 2 1 6-311+G* HF
add_test(NAME HCl_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.75722424")

# test-energy: HCl 0 1 6-311++G* HF
add_test(NAME HCl_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.08640761")

# test-gradient: HCl 0 1 6-311++G* HF
add_test(NAME HCl_0-1_6-311ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(HCl_0-1_6-311ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.56e-02")

# test-energy: HCl 2 1 6-311++G* HF
add_test(NAME HCl_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.75781870")

# test-energy: HCl 0 1 6-311+G** HF
add_test(NAME HCl_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.09416188")

# test-energy: HCl 2 1 6-311+G** HF
add_test(NAME HCl_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.76556159")

# test-energy: HCl 0 1 6-311++G** HF
add_test(NAME HCl_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.09432952")

# test-energy: HCl 2 1 6-311++G** HF
add_test(NAME HCl_2-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.76599466")

# test-energy: HCl 0 1 CC-PVDZ HF
add_test(NAME HCl_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -460.08865013")

# test-energy: HCl 0 1 CC-PVDZ MP2
add_test(NAME HCl_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(HCl_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -460.24160546")

# test-energy: HCl 1 2 CC-PVDZ HF
add_test(NAME HCl_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.66504132")

# test-energy: HCl -1 2 CC-PVDZ HF
add_test(NAME HCl_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -459.96658886")

# test-energy: HCl 2 1 CC-PVDZ HF
add_test(NAME HCl_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/HCl.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(HCl_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -458.75938580")

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

# test-energy: water 2 1 6-311++G* HF
add_test(NAME water_2-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.43735855")

# test-energy: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04815693")

# test-gradient: water 0 1 6-311+G** HF
add_test(NAME water_0-1_6-311pgss_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-311pgss_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.60e-02")

# test-energy: water 2 1 6-311+G** HF
add_test(NAME water_2-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_2-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.45108756")

# test-energy: water 0 1 6-311++G** HF
add_test(NAME water_0-1_6-311ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-311ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.04827538")

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

# test-energy: formaldehyde 0 1 MINI HF
add_test(NAME formaldehyde_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.97366078")

# test-energy: formaldehyde 0 1 MINI MP2
add_test(NAME formaldehyde_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.07721539")

# test-gradient: formaldehyde 0 1 MINI HF
add_test(NAME formaldehyde_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.03e-01")

# test-energy: formaldehyde 1 2 MINI HF
add_test(NAME formaldehyde_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.56956100")

# test-energy: formaldehyde -1 2 MINI HF
add_test(NAME formaldehyde_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.90666146")

# test-energy: formaldehyde 2 1 MINI HF
add_test(NAME formaldehyde_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -111.65598401")

# test-energy: formaldehyde 0 1 MIDI HF
add_test(NAME formaldehyde_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.17647862")

# test-energy: formaldehyde 0 1 MIDI MP2
add_test(NAME formaldehyde_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.39762378")

# test-gradient: formaldehyde 0 1 MIDI HF
add_test(NAME formaldehyde_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.35e-02")

# test-energy: formaldehyde 1 2 MIDI HF
add_test(NAME formaldehyde_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.83313377")

# test-energy: formaldehyde -1 2 MIDI HF
add_test(NAME formaldehyde_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.07637695")

# test-energy: formaldehyde 2 1 MIDI HF
add_test(NAME formaldehyde_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -111.98330348")

# test-energy: formaldehyde 0 1 STO-2G HF
add_test(NAME formaldehyde_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -109.02412376")

# test-energy: formaldehyde 0 1 STO-2G MP2
add_test(NAME formaldehyde_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -109.13824688")

# test-gradient: formaldehyde 0 1 STO-2G HF
add_test(NAME formaldehyde_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.65e-02")

# test-energy: formaldehyde 1 2 STO-2G HF
add_test(NAME formaldehyde_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -108.79336201")

# test-energy: formaldehyde -1 2 STO-2G HF
add_test(NAME formaldehyde_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -108.72713596")

# test-energy: formaldehyde 2 1 STO-2G HF
add_test(NAME formaldehyde_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -107.98744892")

# test-energy: formaldehyde 0 1 STO-3G HF
add_test(NAME formaldehyde_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.35434714")

# test-energy: formaldehyde 0 1 STO-3G MP2
add_test(NAME formaldehyde_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -112.47215594")

# test-gradient: formaldehyde 0 1 STO-3G HF
add_test(NAME formaldehyde_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.15e-05")

# test-energy: formaldehyde 1 2 STO-3G HF
add_test(NAME formaldehyde_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.08575682")

# test-energy: formaldehyde -1 2 STO-3G HF
add_test(NAME formaldehyde_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.10563806")

# test-energy: formaldehyde 2 1 STO-3G HF
add_test(NAME formaldehyde_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -111.24054424")

# test-energy: formaldehyde 0 1 STO-4G HF
add_test(NAME formaldehyde_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.16104588")

# test-energy: formaldehyde 0 1 STO-4G MP2
add_test(NAME formaldehyde_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.27984156")

# test-gradient: formaldehyde 0 1 STO-4G HF
add_test(NAME formaldehyde_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.51e-03")

# test-energy: formaldehyde 1 2 STO-4G HF
add_test(NAME formaldehyde_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.88836112")

# test-energy: formaldehyde -1 2 STO-4G HF
add_test(NAME formaldehyde_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.91920655")

# test-energy: formaldehyde 2 1 STO-4G HF
add_test(NAME formaldehyde_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.03805728")

# test-energy: formaldehyde 0 1 STO-5G HF
add_test(NAME formaldehyde_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.37512601")

# test-energy: formaldehyde 0 1 STO-5G MP2
add_test(NAME formaldehyde_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.49401502")

# test-gradient: formaldehyde 0 1 STO-5G HF
add_test(NAME formaldehyde_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.25e-03")

# test-energy: formaldehyde 1 2 STO-5G HF
add_test(NAME formaldehyde_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.10110499")

# test-energy: formaldehyde -1 2 STO-5G HF
add_test(NAME formaldehyde_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.13486311")

# test-energy: formaldehyde 2 1 STO-5G HF
add_test(NAME formaldehyde_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.24924757")

# test-energy: formaldehyde 0 1 STO-6G HF
add_test(NAME formaldehyde_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.44075377")

# test-energy: formaldehyde 0 1 STO-6G MP2
add_test(NAME formaldehyde_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.55969222")

# test-gradient: formaldehyde 0 1 STO-6G HF
add_test(NAME formaldehyde_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.85e-03")

# test-energy: formaldehyde 1 2 STO-6G HF
add_test(NAME formaldehyde_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.16650333")

# test-energy: formaldehyde -1 2 STO-6G HF
add_test(NAME formaldehyde_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.20074189")

# test-energy: formaldehyde 2 1 STO-6G HF
add_test(NAME formaldehyde_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.31438858")

# test-energy: formaldehyde 0 1 3-21G HF
add_test(NAME formaldehyde_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.22119680")

# test-energy: formaldehyde 0 1 3-21G MP2
add_test(NAME formaldehyde_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.44180884")

# test-gradient: formaldehyde 0 1 3-21G HF
add_test(NAME formaldehyde_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.68e-02")

# test-energy: formaldehyde 1 2 3-21G HF
add_test(NAME formaldehyde_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.87770124")

# test-energy: formaldehyde -1 2 3-21G HF
add_test(NAME formaldehyde_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.11933362")

# test-energy: formaldehyde 2 1 3-21G HF
add_test(NAME formaldehyde_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.02639053")

# test-energy: formaldehyde 0 1 4-31G HF
add_test(NAME formaldehyde_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.69171435")

# test-energy: formaldehyde 0 1 4-31G MP2
add_test(NAME formaldehyde_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.91523506")

# test-gradient: formaldehyde 0 1 4-31G HF
add_test(NAME formaldehyde_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.08e-02")

# test-energy: formaldehyde 1 2 4-31G HF
add_test(NAME formaldehyde_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.34489975")

# test-energy: formaldehyde -1 2 4-31G HF
add_test(NAME formaldehyde_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.60635179")

# test-energy: formaldehyde 2 1 4-31G HF
add_test(NAME formaldehyde_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.49499039")

# test-energy: formaldehyde 0 1 6-21G HF
add_test(NAME formaldehyde_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.69799625")

# test-energy: formaldehyde 0 1 6-21G MP2
add_test(NAME formaldehyde_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.91828761")

# test-gradient: formaldehyde 0 1 6-21G HF
add_test(NAME formaldehyde_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.34e-02")

# test-energy: formaldehyde 1 2 6-21G HF
add_test(NAME formaldehyde_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.35654118")

# test-energy: formaldehyde -1 2 6-21G HF
add_test(NAME formaldehyde_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.59428475")

# test-energy: formaldehyde 2 1 6-21G HF
add_test(NAME formaldehyde_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.50751196")

# test-energy: formaldehyde 0 1 6-31G HF
add_test(NAME formaldehyde_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.80762059")

# test-energy: formaldehyde 0 1 6-31G MP2
add_test(NAME formaldehyde_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -114.02883703")

# test-gradient: formaldehyde 0 1 6-31G HF
add_test(NAME formaldehyde_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.56e-02")

# test-energy: formaldehyde 1 2 6-31G HF
add_test(NAME formaldehyde_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.46024180")

# test-energy: formaldehyde -1 2 6-31G HF
add_test(NAME formaldehyde_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.72660108")

# test-energy: formaldehyde 2 1 6-31G HF
add_test(NAME formaldehyde_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.61147702")

# test-energy: formaldehyde 0 1 6-31+G HF
add_test(NAME formaldehyde_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.81187055")

# test-energy: formaldehyde 0 1 6-31+G MP2
add_test(NAME formaldehyde_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -114.03821974")

# test-gradient: formaldehyde 0 1 6-31+G HF
add_test(NAME formaldehyde_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.59e-02")

# test-energy: formaldehyde 2 1 6-31+G HF
add_test(NAME formaldehyde_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.61485671")

# test-energy: formaldehyde 0 1 6-31++G HF
add_test(NAME formaldehyde_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.81213361")

# test-gradient: formaldehyde 0 1 6-31++G HF
add_test(NAME formaldehyde_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.58e-02")

# test-energy: formaldehyde 2 1 6-31++G HF
add_test(NAME formaldehyde_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.61501070")

# test-energy: formaldehyde 0 1 6-31+G* HF
add_test(NAME formaldehyde_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.86858703")

# test-energy: formaldehyde 2 1 6-31+G* HF
add_test(NAME formaldehyde_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.67149974")

# test-energy: formaldehyde 0 1 6-31++G* HF
add_test(NAME formaldehyde_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.86875341")

# test-energy: formaldehyde 2 1 6-31++G* HF
add_test(NAME formaldehyde_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.67167215")

# test-energy: formaldehyde 0 1 6-311G HF
add_test(NAME formaldehyde_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.83869423")

# test-energy: formaldehyde 2 1 6-311G HF
add_test(NAME formaldehyde_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.64313320")

# test-energy: formaldehyde 2 1 6-311+G HF
add_test(NAME formaldehyde_2-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.64544643")

# test-energy: formaldehyde 0 1 CC-PVDZ HF
add_test(NAME formaldehyde_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.87495844")

# test-energy: formaldehyde 2 1 CC-PVDZ HF
add_test(NAME formaldehyde_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.68502515")

# test-energy: methane 0 1 MINI HF
add_test(NAME methane_0-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.81350709")

# test-energy: methane 0 1 MINI MP2
add_test(NAME methane_0-1_mini_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_mini_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -39.86234767")

# test-gradient: methane 0 1 MINI HF
add_test(NAME methane_0-1_mini_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_mini_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.54e-01")

# test-energy: methane 1 2 MINI HF
add_test(NAME methane_1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.24302783")

# test-energy: methane -1 2 MINI HF
add_test(NAME methane_-1-2_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.36486561")

# test-energy: methane 2 1 MINI HF
add_test(NAME methane_2-1_mini_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "mini" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_mini_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.23810751")

# test-energy: methane 0 1 MIDI HF
add_test(NAME methane_0-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.95639725")

# test-energy: methane 0 1 MIDI MP2
add_test(NAME methane_0-1_midi_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_midi_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.05567338")

# test-gradient: methane 0 1 MIDI HF
add_test(NAME methane_0-1_midi_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_midi_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.81e-03")

# test-energy: methane 1 2 MIDI HF
add_test(NAME methane_1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.45941640")

# test-energy: methane -1 2 MIDI HF
add_test(NAME methane_-1-2_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.72344160")

# test-energy: methane 2 1 MIDI HF
add_test(NAME methane_2-1_midi_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "midi" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_midi_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.49117263")

# test-energy: methane 0 1 STO-2G HF
add_test(NAME methane_0-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.58821853")

# test-energy: methane 0 1 STO-2G MP2
add_test(NAME methane_0-1_sto-2g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-2g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -38.64085000")

# test-gradient: methane 0 1 STO-2G HF
add_test(NAME methane_0-1_sto-2g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-2g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.15e-02")

# test-energy: methane 1 2 STO-2G HF
add_test(NAME methane_1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.10586557")

# test-energy: methane -1 2 STO-2G HF
add_test(NAME methane_-1-2_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -37.86375976")

# test-energy: methane 2 1 STO-2G HF
add_test(NAME methane_2-1_sto-2g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-2g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_sto-2g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -37.13056874")

# test-energy: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.72686367")

# test-energy: methane 0 1 STO-3G MP2
add_test(NAME methane_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -39.78270593")

# test-gradient: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-3g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.23e-06")

# test-energy: methane 1 2 STO-3G HF
add_test(NAME methane_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.23858013")

# test-energy: methane -1 2 STO-3G HF
add_test(NAME methane_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.04630900")

# test-energy: methane 2 1 STO-3G HF
add_test(NAME methane_2-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.25601989")

# test-energy: methane 0 1 STO-4G HF
add_test(NAME methane_0-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.00773213")

# test-energy: methane 0 1 STO-4G MP2
add_test(NAME methane_0-1_sto-4g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-4g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.06437161")

# test-gradient: methane 0 1 STO-4G HF
add_test(NAME methane_0-1_sto-4g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-4g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.83e-03")

# test-energy: methane 1 2 STO-4G HF
add_test(NAME methane_1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.51847592")

# test-energy: methane -1 2 STO-4G HF
add_test(NAME methane_-1-2_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.33327117")

# test-energy: methane 2 1 STO-4G HF
add_test(NAME methane_2-1_sto-4g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-4g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_sto-4g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.53439488")

# test-energy: methane 0 1 STO-5G HF
add_test(NAME methane_0-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.08640356")

# test-energy: methane 0 1 STO-5G MP2
add_test(NAME methane_0-1_sto-5g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-5g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.14310682")

# test-gradient: methane 0 1 STO-5G HF
add_test(NAME methane_0-1_sto-5g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-5g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.87e-03")

# test-energy: methane 1 2 STO-5G HF
add_test(NAME methane_1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.59581607")

# test-energy: methane -1 2 STO-5G HF
add_test(NAME methane_-1-2_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.41326712")

# test-energy: methane 2 1 STO-5G HF
add_test(NAME methane_2-1_sto-5g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-5g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_sto-5g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.61043850")

# test-energy: methane 0 1 STO-6G HF
add_test(NAME methane_0-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.11061063")

# test-energy: methane 0 1 STO-6G MP2
add_test(NAME methane_0-1_sto-6g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-6g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.16734630")

# test-gradient: methane 0 1 STO-6G HF
add_test(NAME methane_0-1_sto-6g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-6g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.39e-03")

# test-energy: methane 1 2 STO-6G HF
add_test(NAME methane_1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.61973654")

# test-energy: methane -1 2 STO-6G HF
add_test(NAME methane_-1-2_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.43821026")

# test-energy: methane 2 1 STO-6G HF
add_test(NAME methane_2-1_sto-6g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-6g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_sto-6g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.63408702")

# test-energy: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.97687752")

# test-energy: methane 0 1 3-21G MP2
add_test(NAME methane_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.07687501")

# test-gradient: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_3-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.59e-04")

# test-energy: methane 1 2 3-21G HF
add_test(NAME methane_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.48044515")

# test-energy: methane -1 2 3-21G HF
add_test(NAME methane_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.69076294")

# test-energy: methane 2 1 3-21G HF
add_test(NAME methane_2-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.50887351")

# test-energy: methane 0 1 4-31G HF
add_test(NAME methane_0-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.13975640")

# test-energy: methane 0 1 4-31G MP2
add_test(NAME methane_0-1_4-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_4-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.24020998")

# test-gradient: methane 0 1 4-31G HF
add_test(NAME methane_0-1_4-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_4-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.81e-03")

# test-energy: methane 1 2 4-31G HF
add_test(NAME methane_1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.64552638")

# test-energy: methane -1 2 4-31G HF
add_test(NAME methane_-1-2_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.88991448")

# test-energy: methane 2 1 4-31G HF
add_test(NAME methane_2-1_4-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "4-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_4-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.67731412")

# test-energy: methane 0 1 6-21G HF
add_test(NAME methane_0-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.14926306")

# test-energy: methane 0 1 6-21G MP2
add_test(NAME methane_0-1_6-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.24914745")

# test-gradient: methane 0 1 6-21G HF
add_test(NAME methane_0-1_6-21g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-21g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.33e-04")

# test-energy: methane 1 2 6-21G HF
add_test(NAME methane_1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.65177572")

# test-energy: methane -1 2 6-21G HF
add_test(NAME methane_-1-2_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.86452249")

# test-energy: methane 2 1 6-21G HF
add_test(NAME methane_2-1_6-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-21g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.68331594")

# test-energy: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18055192")

# test-energy: methane 0 1 6-31G MP2
add_test(NAME methane_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.28000276")

# test-gradient: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.32e-03")

# test-energy: methane 1 2 6-31G HF
add_test(NAME methane_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.68642404")

# test-energy: methane -1 2 6-31G HF
add_test(NAME methane_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.93284040")

# test-energy: methane 2 1 6-31G HF
add_test(NAME methane_2-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.71876086")

# test-energy: methane 0 1 6-31+G HF
add_test(NAME methane_0-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18119106")

# test-energy: methane 0 1 6-31+G MP2
add_test(NAME methane_0-1_6-31pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.28164604")

# test-gradient: methane 0 1 6-31+G HF
add_test(NAME methane_0-1_6-31pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.89e-04")

# test-energy: methane 1 2 6-31+G HF
add_test(NAME methane_1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.68669391")

# test-energy: methane -1 2 6-31+G HF
add_test(NAME methane_-1-2_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.08942753")

# test-energy: methane 2 1 6-31+G HF
add_test(NAME methane_2-1_6-31pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.71927519")

# test-energy: methane 0 1 6-31++G HF
add_test(NAME methane_0-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18121968")

# test-energy: methane 0 1 6-31++G MP2
add_test(NAME methane_0-1_6-31ppg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31ppg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.28189720")

# test-gradient: methane 0 1 6-31++G HF
add_test(NAME methane_0-1_6-31ppg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31ppg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.35e-04")

# test-energy: methane -1 2 6-31++G HF
add_test(NAME methane_-1-2_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.13788454")

# test-energy: methane 2 1 6-31++G HF
add_test(NAME methane_2-1_6-31ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.71957358")

# test-energy: methane 0 1 6-31+G* HF
add_test(NAME methane_0-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.19534781")

# test-energy: methane 0 1 6-31+G* MP2
add_test(NAME methane_0-1_6-31pgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31pgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.33487991")

# test-gradient: methane 0 1 6-31+G* HF
add_test(NAME methane_0-1_6-31pgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31pgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.50e-03")

# test-energy: methane 2 1 6-31+G* HF
add_test(NAME methane_2-1_6-31pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.75744947")

# test-energy: methane 0 1 6-31++G* HF
add_test(NAME methane_0-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.19542981")

# test-energy: methane 0 1 6-31++G* MP2
add_test(NAME methane_0-1_6-31ppgs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31ppgs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.33515399")

# test-gradient: methane 0 1 6-31++G* HF
add_test(NAME methane_0-1_6-31ppgs_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31ppgs_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.56e-03")

# test-energy: methane 2 1 6-31++G* HF
add_test(NAME methane_2-1_6-31ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.75780049")

# test-energy: methane 0 1 6-31+G** HF
add_test(NAME methane_0-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.20210226")

# test-energy: methane 2 1 6-31+G** HF
add_test(NAME methane_2-1_6-31pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31+g**" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-31pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.76848215")

# test-energy: methane 0 1 6-31++G** HF
add_test(NAME methane_0-1_6-31ppgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31++g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31ppgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.20214026")

# test-energy: methane 0 1 6-311G HF
add_test(NAME methane_0-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18817805")

# test-energy: methane 0 1 6-311G MP2
add_test(NAME methane_0-1_6-311g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-311g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.30900220")

# test-gradient: methane 0 1 6-311G HF
add_test(NAME methane_0-1_6-311g_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-311g_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.63e-03")

# test-energy: methane -1 2 6-311G HF
add_test(NAME methane_-1-2_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.03073177")

# test-energy: methane 2 1 6-311G HF
add_test(NAME methane_2-1_6-311g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-311g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.72216325")

# test-energy: methane 0 1 6-311+G HF
add_test(NAME methane_0-1_6-311pg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311pg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18846697")

# test-energy: methane 0 1 6-311+G MP2
add_test(NAME methane_0-1_6-311pg_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-311pg_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.30965001")

# test-gradient: methane 0 1 6-311+G HF
add_test(NAME methane_0-1_6-311pg_hf_gradient COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-311pg_hf_gradient PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.15e-03")

# test-energy: methane 0 1 6-311++G HF
add_test(NAME methane_0-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18850086")

# test-energy: methane 2 1 6-311++G HF
add_test(NAME methane_2-1_6-311ppg_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-311ppg_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.72277575")

# test-energy: methane 0 1 6-311+G* HF
add_test(NAME methane_0-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.20276386")

# test-energy: methane 2 1 6-311+G* HF
add_test(NAME methane_2-1_6-311pgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g*" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_6-311pgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.73413470")

# test-energy: methane 0 1 6-311++G* HF
add_test(NAME methane_0-1_6-311ppgs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311++g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311ppgs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.20290954")

# test-energy: methane 0 1 6-311+G** HF
add_test(NAME methane_0-1_6-311pgss_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-311+g**" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-311pgss_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.20908236")

# test-energy: methane 0 1 CC-PVDZ HF
add_test(NAME methane_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.19854469")

# test-energy: methane 2 1 CC-PVDZ HF
add_test(NAME methane_2-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 2 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_2-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -38.76711423")
