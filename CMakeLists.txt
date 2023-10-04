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
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -static -static-libstdc++ -static-libgcc")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -static -static-libstdc++ -static-libgcc")
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

# find packages
find_package(OpenMP REQUIRED)

# add hazel executable
add_executable(hazel src/ci.cpp src/distributor.cpp src/eigen.cpp src/gradient.cpp src/hessian.cpp src/hf.cpp src/integral.cpp src/lambda.cpp src/main.cpp src/md.cpp src/mp.cpp src/optimizer.cpp src/qd.cpp src/system.cpp src/timer.cpp src/transform.cpp src/utility.cpp)
target_compile_definitions(hazel PRIVATE CXXFLAGS="${CXXFLAGS}" DATADIR="${PROJECT_SOURCE_DIR}")
target_link_libraries(hazel Eigen3::Eigen)

# link libint
if (${STANDALONE})
    target_link_libraries(hazel libint2)
else()
    target_link_libraries(hazel int2.a)
endif()


# link OpenMP on build type
if (CMAKE_BUILD_TYPE STREQUAL Release AND OPENMP_FOUND)
    target_link_libraries(hazel OpenMP::OpenMP_CXX)
endif()

# enable testing
enable_testing()

# add build target before all tests
add_test(build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR})

# TESTS ==============================================================================================================================================

# test-energy: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.58368698")

# test-energy: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.96590121")

# test-energy: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.97974730")

# test-energy: water 0 1 6-31G* HF
add_test(NAME water_0-1_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.00543857")

# test-energy: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(water_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -76.02312305")

# test-energy: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19927942")

# test-energy: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.66228146")

# test-energy: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57902812")

# test-energy: water 1 2 6-31G* HF
add_test(NAME water_1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.61027954")

# test-energy: water 1 2 CC-PVDZ HF
add_test(NAME water_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.63182597")

# test-energy: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.34401430")

# test-energy: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.30165233")

# test-energy: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.79551492")

# test-energy: water -1 2 6-31G* HF
add_test(NAME water_-1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.81088399")

# test-energy: water -1 2 CC-PVDZ HF
add_test(NAME water_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.85263322")

# test-energy: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.86784138")

# test-energy: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.45541979")

# test-energy: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -56.15690488")

# test-energy: ammonia 0 1 6-31G* HF
add_test(NAME ammonia_0-1_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -56.18130480")

# test-energy: ammonia 0 1 CC-PVDZ HF
add_test(NAME ammonia_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -56.19421058")

# test-energy: ammonia 1 2 3-21G HF
add_test(NAME ammonia_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.53488341")

# test-energy: ammonia 1 2 STO-3G HF
add_test(NAME ammonia_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.16926782")

# test-energy: ammonia 1 2 6-31G HF
add_test(NAME ammonia_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.81530957")

# test-energy: ammonia 1 2 6-31G* HF
add_test(NAME ammonia_1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.83383906")

# test-energy: ammonia 1 2 CC-PVDZ HF
add_test(NAME ammonia_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.85009309")

# test-energy: ammonia -1 2 3-21G HF
add_test(NAME ammonia_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.61977173")

# test-energy: ammonia -1 2 STO-3G HF
add_test(NAME ammonia_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -54.85273604")

# test-energy: ammonia -1 2 6-31G HF
add_test(NAME ammonia_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.95770821")

# test-energy: ammonia -1 2 6-31G* HF
add_test(NAME ammonia_-1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.97422353")

# test-energy: ammonia -1 2 CC-PVDZ HF
add_test(NAME ammonia_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -56.01900567")

# test-energy: formaldehyde 0 1 3-21G HF
add_test(NAME formaldehyde_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.22119680")

# test-energy: formaldehyde 0 1 STO-3G HF
add_test(NAME formaldehyde_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.35434714")

# test-energy: formaldehyde 0 1 6-31G HF
add_test(NAME formaldehyde_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.80762059")

# test-energy: formaldehyde 0 1 6-31G* HF
add_test(NAME formaldehyde_0-1_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.86345372")

# test-energy: formaldehyde 0 1 CC-PVDZ HF
add_test(NAME formaldehyde_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.87495844")

# test-energy: formaldehyde 1 2 3-21G HF
add_test(NAME formaldehyde_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.87770124")

# test-energy: formaldehyde 1 2 STO-3G HF
add_test(NAME formaldehyde_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.08575682")

# test-energy: formaldehyde 1 2 6-31G HF
add_test(NAME formaldehyde_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.46024180")

# test-energy: formaldehyde 1 2 6-31G* HF
add_test(NAME formaldehyde_1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.51789635")

# test-energy: formaldehyde 1 2 CC-PVDZ HF
add_test(NAME formaldehyde_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.53180482")

# test-energy: formaldehyde -1 2 3-21G HF
add_test(NAME formaldehyde_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.11933362")

# test-energy: formaldehyde -1 2 STO-3G HF
add_test(NAME formaldehyde_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -112.10563806")

# test-energy: formaldehyde -1 2 6-31G HF
add_test(NAME formaldehyde_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.72660108")

# test-energy: formaldehyde -1 2 6-31G* HF
add_test(NAME formaldehyde_-1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.77131901")

# test-energy: formaldehyde -1 2 CC-PVDZ HF
add_test(NAME formaldehyde_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(formaldehyde_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -113.78583444")

# test-energy: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.97687752")

# test-energy: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.72686367")

# test-energy: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18055192")

# test-energy: methane 0 1 6-31G* HF
add_test(NAME methane_0-1_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.19483930")

# test-energy: methane 0 1 CC-PVDZ HF
add_test(NAME methane_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.19854469")

# test-energy: methane 1 2 3-21G HF
add_test(NAME methane_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.48044515")

# test-energy: methane 1 2 STO-3G HF
add_test(NAME methane_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.23858013")

# test-energy: methane 1 2 6-31G HF
add_test(NAME methane_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.68642404")

# test-energy: methane 1 2 6-31G* HF
add_test(NAME methane_1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.69956332")

# test-energy: methane 1 2 CC-PVDZ HF
add_test(NAME methane_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.70786463")

# test-energy: methane -1 2 3-21G HF
add_test(NAME methane_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.69076294")

# test-energy: methane -1 2 STO-3G HF
add_test(NAME methane_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.04630900")

# test-energy: methane -1 2 6-31G HF
add_test(NAME methane_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.93284040")

# test-energy: methane -1 2 6-31G* HF
add_test(NAME methane_-1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.94580055")

# test-energy: methane -1 2 CC-PVDZ HF
add_test(NAME methane_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.01026779")

# test-energy: ethane 0 1 3-21G HF
add_test(NAME ethane_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.79391796")

# test-energy: ethane 0 1 STO-3G HF
add_test(NAME ethane_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.30617964")

# test-energy: ethane 0 1 6-31G HF
add_test(NAME ethane_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -79.19748091")

# test-energy: ethane 0 1 6-31G* HF
add_test(NAME ethane_0-1_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_0-1_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -79.22815454")

# test-energy: ethane 0 1 CC-PVDZ HF
add_test(NAME ethane_0-1_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_0-1_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -79.23464790")

# test-energy: ethane 1 2 3-21G HF
add_test(NAME ethane_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.33783478")

# test-energy: ethane 1 2 STO-3G HF
add_test(NAME ethane_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -77.88167031")

# test-energy: ethane 1 2 6-31G HF
add_test(NAME ethane_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.74256982")

# test-energy: ethane 1 2 6-31G* HF
add_test(NAME ethane_1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.79136147")

# test-energy: ethane 1 2 CC-PVDZ HF
add_test(NAME ethane_1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.80043936")

# test-energy: ethane -1 2 3-21G HF
add_test(NAME ethane_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.52357955")

# test-energy: ethane -1 2 STO-3G HF
add_test(NAME ethane_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -77.58461703")

# test-energy: ethane -1 2 6-31G HF
add_test(NAME ethane_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.96385519")

# test-energy: ethane -1 2 6-31G* HF
add_test(NAME ethane_-1-2_6-31gs_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_-1-2_6-31gs_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -78.99305962")

# test-energy: ethane -1 2 CC-PVDZ HF
add_test(NAME ethane_-1-2_cc-pvdz_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ethane_-1-2_cc-pvdz_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -79.05315489")

# test-energy: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.70882320")

# test-energy: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.00485494")

# test-energy: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.11188710")

# test-energy: water 0 1 6-31G* MP2
add_test(NAME water_0-1_6-31gs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31gs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.19454516")

# test-energy: water 0 1 CC-PVDZ MP2
add_test(NAME water_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.22985471")

# test-energy: ammonia 0 1 3-21G MP2
add_test(NAME ammonia_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -55.98296808")

# test-energy: ammonia 0 1 STO-3G MP2
add_test(NAME ammonia_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -55.50545856")

# test-energy: ammonia 0 1 6-31G MP2
add_test(NAME ammonia_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -56.27598095")

# test-energy: ammonia 0 1 6-31G* MP2
add_test(NAME ammonia_0-1_6-31gs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_6-31gs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -56.35320162")

# test-energy: ammonia 0 1 CC-PVDZ MP2
add_test(NAME ammonia_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -56.38485194")

# test-energy: formaldehyde 0 1 3-21G MP2
add_test(NAME formaldehyde_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -113.44180884")

# test-energy: formaldehyde 0 1 STO-3G MP2
add_test(NAME formaldehyde_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -112.47215594")

# test-energy: formaldehyde 0 1 6-31G MP2
add_test(NAME formaldehyde_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -114.02883703")

# test-energy: formaldehyde 0 1 6-31G* MP2
add_test(NAME formaldehyde_0-1_6-31gs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_6-31gs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -114.16765298")

# test-energy: formaldehyde 0 1 CC-PVDZ MP2
add_test(NAME formaldehyde_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(formaldehyde_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -114.19754586")

# test-energy: methane 0 1 3-21G MP2
add_test(NAME methane_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.07687501")

# test-energy: methane 0 1 STO-3G MP2
add_test(NAME methane_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -39.78270593")

# test-energy: methane 0 1 6-31G MP2
add_test(NAME methane_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.28000276")

# test-energy: methane 0 1 6-31G* MP2
add_test(NAME methane_0-1_6-31gs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31gs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.33335608")

# test-energy: methane 0 1 CC-PVDZ MP2
add_test(NAME methane_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.36238873")

# test-energy: ethane 0 1 3-21G MP2
add_test(NAME ethane_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ethane_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -78.98490969")

# test-energy: ethane 0 1 STO-3G MP2
add_test(NAME ethane_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ethane_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -78.41271614")

# test-energy: ethane 0 1 6-31G MP2
add_test(NAME ethane_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ethane_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -79.38752374")

# test-energy: ethane 0 1 6-31G* MP2
add_test(NAME ethane_0-1_6-31gs_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ethane_0-1_6-31gs_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -79.49653761")

# test-energy: ethane 0 1 CC-PVDZ MP2
add_test(NAME ethane_0-1_cc-pvdz_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ethane_0-1_cc-pvdz_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -79.54190007")

# test-energy: water 0 1 3-21G CID
add_test(NAME water_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -75.61399449")

# test-energy: water 0 1 STO-3G CID
add_test(NAME water_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -74.98727979")

# test-energy: water 0 1 6-31G CID
add_test(NAME water_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -76.00904138")

# test-energy: water 0 1 6-31G* CID
add_test(NAME water_0-1_6-31gs_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_6-31gs_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -76.04763081")

# test-energy: water 0 1 CC-PVDZ CID
add_test(NAME water_0-1_cc-pvdz_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_cc-pvdz_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -76.06081556")

# test-energy: ammonia 0 1 3-21G CID
add_test(NAME ammonia_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -55.89382350")

# test-energy: ammonia 0 1 STO-3G CID
add_test(NAME ammonia_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -55.47978199")

# test-energy: ammonia 0 1 6-31G CID
add_test(NAME ammonia_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -56.18529779")

# test-energy: ammonia 0 1 6-31G* CID
add_test(NAME ammonia_0-1_6-31gs_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_6-31gs_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -56.21617777")

# test-energy: ammonia 0 1 CC-PVDZ CID
add_test(NAME ammonia_0-1_cc-pvdz_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_cc-pvdz_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -56.22418393")

# test-energy: formaldehyde 0 1 3-21G CID
add_test(NAME formaldehyde_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(formaldehyde_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -113.29812199")

# test-energy: formaldehyde 0 1 STO-3G CID
add_test(NAME formaldehyde_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(formaldehyde_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -112.48850625")

# test-energy: formaldehyde 0 1 6-31G CID
add_test(NAME formaldehyde_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(formaldehyde_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -113.87465946")

# test-energy: formaldehyde 0 1 6-31G* CID
add_test(NAME formaldehyde_0-1_6-31gs_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(formaldehyde_0-1_6-31gs_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -113.92495474")

# test-energy: formaldehyde 0 1 CC-PVDZ CID
add_test(NAME formaldehyde_0-1_cc-pvdz_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(formaldehyde_0-1_cc-pvdz_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -113.93486750")

# test-energy: methane 0 1 3-21G CID
add_test(NAME methane_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -39.99726353")

# test-energy: methane 0 1 STO-3G CID
add_test(NAME methane_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -39.75564906")

# test-energy: methane 0 1 6-31G CID
add_test(NAME methane_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -40.20078189")

# test-energy: methane 0 1 6-31G* CID
add_test(NAME methane_0-1_6-31gs_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_6-31gs_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -40.21893974")

# test-energy: methane 0 1 CC-PVDZ CID
add_test(NAME methane_0-1_cc-pvdz_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_cc-pvdz_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -40.22128770")

# test-energy: ethane 0 1 3-21G CID
add_test(NAME ethane_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ethane_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -78.81153583")

# test-energy: ethane 0 1 STO-3G CID
add_test(NAME ethane_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ethane_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -78.32970152")

# test-energy: ethane 0 1 6-31G CID
add_test(NAME ethane_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ethane_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -79.21490514")

# test-energy: ethane 0 1 6-31G* CID
add_test(NAME ethane_0-1_6-31gs_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ethane_0-1_6-31gs_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -79.24876213")

# test-energy: ethane 0 1 CC-PVDZ CID
add_test(NAME ethane_0-1_cc-pvdz_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ethane_0-1_cc-pvdz_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -79.25378699")

# test-angrad: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.89e-02")

# test-angrad: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.59e-06")

# test-angrad: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.24e-02")

# test-angrad: water 0 1 6-31G* HF
add_test(NAME water_0-1_6-31gs_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31gs_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.29e-02")

# test-angrad: water 0 1 CC-PVDZ HF
add_test(NAME water_0-1_cc-pvdz_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_cc-pvdz_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.49e-02")

# test-angrad: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.08e-02")

# test-angrad: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.05e-06")

# test-angrad: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.27e-02")

# test-angrad: ammonia 0 1 6-31G* HF
add_test(NAME ammonia_0-1_6-31gs_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_6-31gs_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.42e-02")

# test-angrad: ammonia 0 1 CC-PVDZ HF
add_test(NAME ammonia_0-1_cc-pvdz_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_cc-pvdz_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.43e-02")

# test-angrad: formaldehyde 0 1 3-21G HF
add_test(NAME formaldehyde_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.68e-02")

# test-angrad: formaldehyde 0 1 STO-3G HF
add_test(NAME formaldehyde_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.15e-05")

# test-angrad: formaldehyde 0 1 6-31G HF
add_test(NAME formaldehyde_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.56e-02")

# test-angrad: formaldehyde 0 1 6-31G* HF
add_test(NAME formaldehyde_0-1_6-31gs_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_6-31gs_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.98e-02")

# test-angrad: formaldehyde 0 1 CC-PVDZ HF
add_test(NAME formaldehyde_0-1_cc-pvdz_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/formaldehyde.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(formaldehyde_0-1_cc-pvdz_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.56e-02")

# test-angrad: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.59e-04")

# test-angrad: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.23e-06")

# test-angrad: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.32e-03")

# test-angrad: methane 0 1 6-31G* HF
add_test(NAME methane_0-1_6-31gs_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31gs_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.04e-03")

# test-angrad: methane 0 1 CC-PVDZ HF
add_test(NAME methane_0-1_cc-pvdz_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_cc-pvdz_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.15e-02")

# test-angrad: ethane 0 1 3-21G HF
add_test(NAME ethane_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ethane_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.53e-03")

# test-angrad: ethane 0 1 STO-3G HF
add_test(NAME ethane_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ethane_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.25e-05")

# test-angrad: ethane 0 1 6-31G HF
add_test(NAME ethane_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ethane_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.81e-03")

# test-angrad: ethane 0 1 6-31G* HF
add_test(NAME ethane_0-1_6-31gs_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g*" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ethane_0-1_6-31gs_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.03e-02")

# test-angrad: ethane 0 1 CC-PVDZ HF
add_test(NAME ethane_0-1_cc-pvdz_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "cc-pvdz" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ethane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ethane_0-1_cc-pvdz_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.17e-02")
