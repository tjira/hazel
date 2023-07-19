#include "../include/gradient.h"

int test_grad_ethylene_hf_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "3-21G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.01819013524635, 0.00024787079493, 0.00485083129517, -0.01820703117586, -0.00025627949988, -0.00479000652243, -0.00319069635535, -0.00476019863803, -0.00160339373163, -0.00368241576347, 0.00466966468791, -0.00024970776440, 0.00368937321868, -0.00466522082114, 0.00021892585196, 0.00320063483119, 0.00476416347616, 0.00157335087126;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
