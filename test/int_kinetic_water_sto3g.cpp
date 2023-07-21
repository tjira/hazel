#include "../include/integral.h"

int test_int_kinetic_water_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 29.00320406467810, -0.16801096113783, -0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00454750335691, -0.00454759316824, -0.16801096113783, 0.80812790277365, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.11368846936531, 0.11368775772077, -0.00000000000000, -0.00000000000000, 2.52873122631647, 0.00000000000000, 0.00000000000000, 0.18295676077128, 0.11248713152950, -0.00000000000000, -0.00000000000000, 0.00000000000000, 2.52873122631647, 0.00000000000000, -0.13186292368420, 0.23421773019345, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.52873122631647, 0.12953857981697, -0.01137291156241, -0.00454750335691, 0.11368846936531, 0.18295676077128, -0.13186292368420, 0.12953857981697, 0.76003187992239, 0.00830169969033, -0.00454759316824, 0.11368775772077, 0.11248713152950, 0.23421773019345, -0.01137291156241, 0.00830169969033, 0.76003187992239;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
