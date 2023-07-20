#include "../include/gradient.h"

int test_grad_ethane_hf_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "CC-PVDZ", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << -0.01298130182920, 0.00038805275048, 0.00035999774066, 0.01298129236349, -0.00038869508382, -0.00036064351804, -0.00267587955400, -0.00004113722206, -0.00389560793699, -0.00260526492222, -0.00329975616725, 0.00215941898432, -0.00240606271088, 0.00357137825883, 0.00194995582081, 0.00260533310249, 0.00329910978844, -0.00216039419700, 0.00240567839905, -0.00357139039213, -0.00194869562320, 0.00267620515133, 0.00004243806749, 0.00389596872945;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
