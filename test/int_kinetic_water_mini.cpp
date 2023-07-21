#include "../include/integral.h"

int test_int_kinetic_water_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 29.55420198324335, -0.93007909770272, -0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00065207625116, -0.00065214406659, -0.93007909770272, 0.85977235460205, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.12047483537691, 0.12047424735577, -0.00000000000000, 0.00000000000000, 2.50867322208213, 0.00000000000000, -0.00000000000000, 0.17302669574951, 0.10638199125513, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.50867322208213, -0.00000000000000, -0.12470600091936, 0.22150576858387, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 2.50867322208213, 0.12250781191866, -0.01075565677538, -0.00065207625116, 0.12047483537691, 0.17302669574951, -0.12470600091936, 0.12250781191866, 0.49698730092503, 0.03603017634332, -0.00065214406659, 0.12047424735577, 0.10638199125513, 0.22150576858387, -0.01075565677538, 0.03603017634332, 0.49698730092503;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
