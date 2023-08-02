#include "../include/integral.h"

int test_integral_kinetic_ammonia_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ammonia.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 21.99075331665703, -0.12841636879070, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00726254571116, -0.00726253636194, -0.00726253769366, -0.12841636879070, 0.60699384694403, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.11477492391082, 0.11477499947276, 0.11477498870956, -0.00000000000000, 0.00000000000000, 1.89935812090028, 0.00000000000000, 0.00000000000000, 0.08905953605172, 0.22365987685743, -0.12614665441397, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.89935812090028, 0.00000000000000, 0.19471728172435, -0.09540202964889, 0.13295935384718, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.89935812090028, 0.16401881023944, -0.11671921531559, -0.19788250183381, -0.00726254571116, 0.11477492391082, 0.08905953605172, 0.19471728172435, 0.16401881023944, 0.76003187992239, 0.00076751071982, 0.00076750635679, -0.00726253636194, 0.11477499947276, 0.22365987685743, -0.09540202964889, -0.11671921531559, 0.00076751071982, 0.76003187992239, 0.00076736425263, -0.00726253769366, 0.11477498870956, -0.12614665441397, 0.13295935384718, -0.19788250183381, 0.00076750635679, 0.00076736425263, 0.76003187992239;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}