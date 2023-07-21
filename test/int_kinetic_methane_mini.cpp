#include "../include/integral.h"

int test_int_kinetic_methane_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 16.19074006781007, -0.48698724806222, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00318102023777, -0.00318102700395, -0.00318103469420, -0.00318103796084, -0.48698724806222, 0.44072788591624, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.12022829170512, 0.12022822757469, 0.12022815468570, 0.12022812372407, -0.00000000000000, 0.00000000000000, 1.23814483672461, 0.00000000000000, 0.00000000000000, 0.11847057501185, 0.04595985122141, 0.07494447554930, -0.23937639953144, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.23814483672461, 0.00000000000000, 0.03607876401579, 0.18460715949363, -0.20865992926108, -0.01202739872469, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.23814483672461, 0.20967162647908, -0.15200984368218, -0.10071452828920, 0.04305017664329, -0.00318102023777, 0.12022829170512, 0.11847057501185, 0.03607876401579, 0.20967162647908, 0.49698730092503, 0.01421359062960, 0.01421387943979, 0.01421376678096, -0.00318102700395, 0.12022822757469, 0.04595985122141, 0.18460715949363, -0.15200984368218, 0.01421359062960, 0.49698730092503, 0.01421414727986, 0.01421412480212, -0.00318103469420, 0.12022815468570, 0.07494447554930, -0.20865992926108, -0.10071452828920, 0.01421387943979, 0.01421414727986, 0.49698730092503, 0.01421411748342, -0.00318103796084, 0.12022812372407, -0.23937639953144, -0.01202739872469, 0.04305017664329, 0.01421376678096, 0.01421412480212, 0.01421411748342, 0.49698730092503;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
