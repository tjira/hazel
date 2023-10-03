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

# test-energy: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.19927942")

# test-energy: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.66228146")

# test-energy: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.57902812")

# test-energy: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.34401430")

# test-energy: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -74.30165233")

# test-energy: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(water_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -75.79551492")

# test-energy: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.86784138")

# test-energy: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.45541979")

# test-energy: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -56.15690488")

# test-energy: ammonia 1 2 3-21G HF
add_test(NAME ammonia_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.53488341")

# test-energy: ammonia 1 2 STO-3G HF
add_test(NAME ammonia_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.16926782")

# test-energy: ammonia 1 2 6-31G HF
add_test(NAME ammonia_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.81530957")

# test-energy: ammonia -1 2 3-21G HF
add_test(NAME ammonia_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.61977173")

# test-energy: ammonia -1 2 STO-3G HF
add_test(NAME ammonia_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -54.85273604")

# test-energy: ammonia -1 2 6-31G HF
add_test(NAME ammonia_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(ammonia_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -55.95770821")

# test-energy: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.97687752")

# test-energy: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.72686367")

# test-energy: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8)
set_tests_properties(methane_0-1_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -40.18055192")

# test-energy: methane 1 2 3-21G HF
add_test(NAME methane_1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.48044515")

# test-energy: methane 1 2 STO-3G HF
add_test(NAME methane_1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.23858013")

# test-energy: methane 1 2 6-31G HF
add_test(NAME methane_1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.68642404")

# test-energy: methane -1 2 3-21G HF
add_test(NAME methane_-1-2_3-21g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_3-21g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.69076294")

# test-energy: methane -1 2 STO-3G HF
add_test(NAME methane_-1-2_sto-3g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_sto-3g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.04630900")

# test-energy: methane -1 2 6-31G HF
add_test(NAME methane_-1-2_6-31g_hf_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8)
set_tests_properties(methane_-1-2_6-31g_hf_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL HARTREE-FOCK ENERGY: -39.93284040")

# test-energy: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.70882320")

# test-energy: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -75.00485494")

# test-energy: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(water_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -76.11188710")

# test-energy: ammonia 0 1 3-21G MP2
add_test(NAME ammonia_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -55.98296808")

# test-energy: ammonia 0 1 STO-3G MP2
add_test(NAME ammonia_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -55.50545856")

# test-energy: ammonia 0 1 6-31G MP2
add_test(NAME ammonia_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(ammonia_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -56.27598095")

# test-energy: methane 0 1 3-21G MP2
add_test(NAME methane_0-1_3-21g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_3-21g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.07687501")

# test-energy: methane 0 1 STO-3G MP2
add_test(NAME methane_0-1_sto-3g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_sto-3g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -39.78270593")

# test-energy: methane 0 1 6-31G MP2
add_test(NAME methane_0-1_6-31g_mp2_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2)
set_tests_properties(methane_0-1_6-31g_mp2_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL MP2 ENERGY: -40.28000276")

# test-energy: water 0 1 3-21G CID
add_test(NAME water_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -75.61399449")

# test-energy: water 0 1 STO-3G CID
add_test(NAME water_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -74.98727979")

# test-energy: water 0 1 6-31G CID
add_test(NAME water_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(water_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -76.00904138")

# test-energy: ammonia 0 1 3-21G CID
add_test(NAME ammonia_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -55.89382350")

# test-energy: ammonia 0 1 STO-3G CID
add_test(NAME ammonia_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -55.47978199")

# test-energy: ammonia 0 1 6-31G CID
add_test(NAME ammonia_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(ammonia_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -56.18529779")

# test-energy: methane 0 1 3-21G CID
add_test(NAME methane_0-1_3-21g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_3-21g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -39.99726353")

# test-energy: methane 0 1 STO-3G CID
add_test(NAME methane_0-1_sto-3g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_sto-3g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -39.75564906")

# test-energy: methane 0 1 6-31G CID
add_test(NAME methane_0-1_6-31g_cid_energy COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 ci -e d)
set_tests_properties(methane_0-1_6-31g_cid_energy PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "FINAL CI ENERGY: -40.20078189")

# test-angrad: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.89e-02")

# test-angrad: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.59e-06")

# test-angrad: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(water_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.24e-02")

# test-angrad: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.08e-02")

# test-angrad: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.05e-06")

# test-angrad: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(ammonia_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.27e-02")

# test-angrad: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_3-21g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.59e-04")

# test-angrad: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_sto-3g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.23e-06")

# test-angrad: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_angrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 0 1e-5)
set_tests_properties(methane_0-1_6-31g_hf_angrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.32e-03")

# test-numgrad: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_0-1_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.88e-02")

# test-numgrad: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_0-1_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.50e-05")

# test-numgrad: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_0-1_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 8.24e-02")

# test-numgrad: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.49e-02")

# test-numgrad: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.03e-02")

# test-numgrad: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.80e-02")

# test-numgrad: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_-1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.43e-01")

# test-numgrad: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_-1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.73e-01")

# test-numgrad: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(water_-1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.32e-02")

# test-numgrad: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_0-1_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.08e-02")

# test-numgrad: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_0-1_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.97e-05")

# test-numgrad: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_0-1_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.26e-02")

# test-numgrad: ammonia 1 2 3-21G HF
add_test(NAME ammonia_1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.35e-01")

# test-numgrad: ammonia 1 2 STO-3G HF
add_test(NAME ammonia_1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.10e-01")

# test-numgrad: ammonia 1 2 6-31G HF
add_test(NAME ammonia_1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.50e-01")

# test-numgrad: ammonia -1 2 3-21G HF
add_test(NAME ammonia_-1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_-1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 9.61e-02")

# test-numgrad: ammonia -1 2 STO-3G HF
add_test(NAME ammonia_-1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_-1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.93e-01")

# test-numgrad: ammonia -1 2 6-31G HF
add_test(NAME ammonia_-1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(ammonia_-1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.84e-02")

# test-numgrad: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_0-1_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.43e-04")

# test-numgrad: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_0-1_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.24e-04")

# test-numgrad: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_0-1_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.24e-03")

# test-numgrad: methane 1 2 3-21G HF
add_test(NAME methane_1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.29e-01")

# test-numgrad: methane 1 2 STO-3G HF
add_test(NAME methane_1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.66e-01")

# test-numgrad: methane 1 2 6-31G HF
add_test(NAME methane_1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.27e-01")

# test-numgrad: methane -1 2 3-21G HF
add_test(NAME methane_-1-2_3-21g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_-1-2_3-21g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.80e-02")

# test-numgrad: methane -1 2 STO-3G HF
add_test(NAME methane_-1-2_sto-3g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_-1-2_sto-3g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 7.40e-01")

# test-numgrad: methane -1 2 6-31G HF
add_test(NAME methane_-1-2_6-31g_hf_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -g 1 1e-5)
set_tests_properties(methane_-1-2_6-31g_hf_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 6.01e-02")

# test-numgrad: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(water_0-1_3-21g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.74e-02")

# test-numgrad: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(water_0-1_sto-3g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 5.13e-02")

# test-numgrad: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(water_0-1_6-31g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 3.68e-02")

# test-numgrad: ammonia 0 1 3-21G MP2
add_test(NAME ammonia_0-1_3-21g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(ammonia_0-1_3-21g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.52e-02")

# test-numgrad: ammonia 0 1 STO-3G MP2
add_test(NAME ammonia_0-1_sto-3g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(ammonia_0-1_sto-3g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.61e-02")

# test-numgrad: ammonia 0 1 6-31G MP2
add_test(NAME ammonia_0-1_6-31g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(ammonia_0-1_6-31g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 4.07e-02")

# test-numgrad: methane 0 1 3-21G MP2
add_test(NAME methane_0-1_3-21g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(methane_0-1_3-21g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.29e-02")

# test-numgrad: methane 0 1 STO-3G MP2
add_test(NAME methane_0-1_sto-3g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(methane_0-1_sto-3g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 2.72e-02")

# test-numgrad: methane 0 1 6-31G MP2
add_test(NAME methane_0-1_6-31g_mp2_numgrad COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -g 1 1e-5)
set_tests_properties(methane_0-1_6-31g_mp2_numgrad PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "GRADIENT NORM: 1.73e-02")

# test-numhess: water 0 1 3-21G HF
add_test(NAME water_0-1_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_0-1_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.34e\\+00")

# test-numhess: water 0 1 STO-3G HF
add_test(NAME water_0-1_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_0-1_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.90e\\+00")

# test-numhess: water 0 1 6-31G HF
add_test(NAME water_0-1_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_0-1_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.29e\\+00")

# test-numhess: water 1 2 3-21G HF
add_test(NAME water_1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.33e\\+00")

# test-numhess: water 1 2 STO-3G HF
add_test(NAME water_1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.96e\\+00")

# test-numhess: water 1 2 6-31G HF
add_test(NAME water_1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.35e\\+00")

# test-numhess: water -1 2 3-21G HF
add_test(NAME water_-1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_-1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 7.86e-01")

# test-numhess: water -1 2 STO-3G HF
add_test(NAME water_-1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_-1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 6.55e\\+00")

# test-numhess: water -1 2 6-31G HF
add_test(NAME water_-1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(water_-1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 7.00e-01")

# test-numhess: ammonia 0 1 3-21G HF
add_test(NAME ammonia_0-1_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_0-1_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.50e\\+00")

# test-numhess: ammonia 0 1 STO-3G HF
add_test(NAME ammonia_0-1_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_0-1_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.12e\\+00")

# test-numhess: ammonia 0 1 6-31G HF
add_test(NAME ammonia_0-1_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_0-1_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.54e\\+00")

# test-numhess: ammonia 1 2 3-21G HF
add_test(NAME ammonia_1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.50e\\+00")

# test-numhess: ammonia 1 2 STO-3G HF
add_test(NAME ammonia_1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.20e\\+00")

# test-numhess: ammonia 1 2 6-31G HF
add_test(NAME ammonia_1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.58e\\+00")

# test-numhess: ammonia -1 2 3-21G HF
add_test(NAME ammonia_-1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_-1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 8.09e-01")

# test-numhess: ammonia -1 2 STO-3G HF
add_test(NAME ammonia_-1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_-1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.18e\\+01")

# test-numhess: ammonia -1 2 6-31G HF
add_test(NAME ammonia_-1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(ammonia_-1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.03e\\+00")

# test-numhess: methane 0 1 3-21G HF
add_test(NAME methane_0-1_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_0-1_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.65e\\+00")

# test-numhess: methane 0 1 STO-3G HF
add_test(NAME methane_0-1_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_0-1_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.24e\\+00")

# test-numhess: methane 0 1 6-31G HF
add_test(NAME methane_0-1_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_0-1_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.65e\\+00")

# test-numhess: methane 1 2 3-21G HF
add_test(NAME methane_1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.37e\\+00")

# test-numhess: methane 1 2 STO-3G HF
add_test(NAME methane_1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 4.33e\\+00")

# test-numhess: methane 1 2 6-31G HF
add_test(NAME methane_1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.45e\\+00")

# test-numhess: methane -1 2 3-21G HF
add_test(NAME methane_-1-2_3-21g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_-1-2_3-21g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 7.09e-01")

# test-numhess: methane -1 2 STO-3G HF
add_test(NAME methane_-1-2_sto-3g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_-1-2_sto-3g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.24e\\+01")

# test-numhess: methane -1 2 6-31G HF
add_test(NAME methane_-1-2_6-31g_hf_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c -1 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 2 hf -i 1000 -t 1e-8 -f 1 1e-4)
set_tests_properties(methane_-1-2_6-31g_hf_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.11e\\+00")

# test-numhess: water 0 1 3-21G MP2
add_test(NAME water_0-1_3-21g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(water_0-1_3-21g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.30e\\+00")

# test-numhess: water 0 1 STO-3G MP2
add_test(NAME water_0-1_sto-3g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(water_0-1_sto-3g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.83e\\+00")

# test-numhess: water 0 1 6-31G MP2
add_test(NAME water_0-1_6-31g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/water.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(water_0-1_6-31g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.28e\\+00")

# test-numhess: ammonia 0 1 3-21G MP2
add_test(NAME ammonia_0-1_3-21g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(ammonia_0-1_3-21g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.46e\\+00")

# test-numhess: ammonia 0 1 STO-3G MP2
add_test(NAME ammonia_0-1_sto-3g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(ammonia_0-1_sto-3g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.04e\\+00")

# test-numhess: ammonia 0 1 6-31G MP2
add_test(NAME ammonia_0-1_6-31g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/ammonia.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(ammonia_0-1_6-31g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.52e\\+00")

# test-numhess: methane 0 1 3-21G MP2
add_test(NAME methane_0-1_3-21g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "3-21g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(methane_0-1_3-21g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.63e\\+00")

# test-numhess: methane 0 1 STO-3G MP2
add_test(NAME methane_0-1_sto-3g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "sto-3g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(methane_0-1_sto-3g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 2.14e\\+00")

# test-numhess: methane 0 1 6-31G MP2
add_test(NAME methane_0-1_6-31g_mp2_numhess COMMAND ${PROJECT_SOURCE_DIR}/bin/hazel -b "6-31g" -c 0 -f ${PROJECT_SOURCE_DIR}/example/molecule/methane.xyz -n 2 -s 1 hf -i 1000 -t 1e-8 mp2 -f 1 1e-4)
set_tests_properties(methane_0-1_6-31g_mp2_numhess PROPERTIES DEPENDS build PASS_REGULAR_EXPRESSION "HESSIAN NORM: 1.64e\\+00")
