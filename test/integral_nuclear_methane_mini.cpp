#include "../include/integral.h"

int test_integral_nuclear_methane_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -35.98926363198078, -3.37397332143948, 0.00000008451237, 0.00000007970951, 0.00000014705967, -1.38895154850919, -1.38895096677571, -1.38895029310884, -1.38894999955263, -3.37397332143948, -6.67224652360179, 0.00000110127662, 0.00000102624863, 0.00000185929926, -3.33882803812536, -3.33882794016884, -3.33882729144039, -3.33882669743688, 0.00000008451237, 0.00000110127662, -6.50167781975504, 0.00000011325658, -0.00000078249079, -1.25178711722716, -0.48562195334986, -0.79188066365041, 2.52930874503558, 0.00000007970951, 0.00000102624863, 0.00000011325658, -6.50167731052889, 0.00000022382676, -0.38121600852203, -1.95060211443628, 2.20475090525815, 0.12708479439048, 0.00000014705967, 0.00000185929926, -0.00000078249079, 0.00000022382676, -6.50167903598022, -2.21543833758957, 1.60617275077839, 1.06417455816139, -0.45487751197563, -1.38895154850919, -3.33882803812536, -1.25178711722716, -0.38121600852203, -2.21543833758957, -4.67395464781818, -1.43147767923098, -1.43148401611713, -1.43148144161980, -1.38895096677571, -3.33882794016884, -0.48562195334986, -1.95060211443628, 1.60617275077839, -1.43147767923098, -4.67395566610659, -1.43149005928172, -1.43148945970380, -1.38895029310884, -3.33882729144039, -0.79188066365041, 2.20475090525815, 1.06417455816139, -1.43148401611713, -1.43149005928172, -4.67395578995577, -1.43148931151880, -1.38894999955263, -3.33882669743688, 2.52930874503558, 0.12708479439048, -0.45487751197563, -1.43148144161980, -1.43148945970380, -1.43148931151880, -4.67395523014843;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}