#include "../include/gradient.h"

int test_gradient_hf_ethylene_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "MINI", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // calculate HF energy and gradient
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    Matrix G = Gradient({}).get(system, rhfres, false);
    libint2::finalize();

    // create the expectation gradient
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.02623301088050, 0.00036299197094, 0.00694641059786, -0.02623234386739, -0.00036295331862, -0.00695041385670, 0.03626700289656, 0.06925390458596, 0.02041880184125, 0.04340009910299, -0.06814744824546, 0.00068585669175, -0.04339977047468, 0.06814745648944, -0.00068368080697, -0.03626799853798, -0.06925395148226, -0.02041697446721;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
