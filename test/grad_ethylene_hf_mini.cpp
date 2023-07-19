#include "../include/gradient.h"

int test_grad_ethylene_hf_mini(int, char**) {
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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.02623301087940, 0.00036299197068, 0.00694641059751, -0.02623234386978, -0.00036295331880, -0.00695041385935, 0.03626700289653, 0.06925390458595, 0.02041880184122, 0.04340009910297, -0.06814744824546, 0.00068585669171, -0.04339977047471, 0.06814745648944, -0.00068368080703, -0.03626799853802, -0.06925395148226, -0.02041697446727;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
