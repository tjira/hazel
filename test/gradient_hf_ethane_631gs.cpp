#include "../include/gradient.h"

int test_gradient_hf_ethane_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "6-31G*", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << -0.00690777624283, 0.00020640064664, 0.00019134872292, 0.00690776388224, -0.00020702631295, -0.00019196277061, -0.00078870466876, 0.00005373444726, 0.00096497414053, -0.00080566210302, 0.00082569192976, -0.00047455361426, -0.00085348650540, -0.00080565975869, -0.00042216739858, 0.00080569264857, -0.00082696376575, 0.00047257549713, 0.00085310234881, 0.00080507075961, 0.00042445161234, 0.00078907064031, -0.00005124794586, -0.00096466618947;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
