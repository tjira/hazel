#include "../include/gradient.h"

int test_gradient_hf_ethylene_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "6-31G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.02798173157854, 0.00038367453866, 0.00744479691665, -0.02799879129680, -0.00039214202204, -0.00738331497262, -0.00277881355683, -0.00471290456436, -0.00148793227537, -0.00326616819106, 0.00463389844814, -0.00014616106594, 0.00327326466556, -0.00462942466655, 0.00011506577319, 0.00278877680022, 0.00471689826615, 0.00145754562420;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
