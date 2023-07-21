#include "../include/integral.h"

int test_integral_overlap_water_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.20304660416675, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.05195118411909, 0.05195102613984, 0.20304660416675, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50621341445783, 0.50621249447803, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.27690885240420, 0.17025223381536, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.19957727016003, 0.35449469839259, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.19605932750433, -0.01721320094270, 0.05195118411909, 0.50621341445783, 0.27690885240420, -0.19957727016003, 0.19605932750433, 1.00000000000000, 0.37104400878062, 0.05195102613984, 0.50621249447803, 0.17025223381536, 0.35449469839259, -0.01721320094270, 0.37104400878062, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}
