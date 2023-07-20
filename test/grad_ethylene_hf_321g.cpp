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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.01819013524436, 0.00024787079516, 0.00485083129513, -0.01820703117555, -0.00025627950003, -0.00479000652241, -0.00319069635536, -0.00476019863804, -0.00160339373162, -0.00368241576349, 0.00466966468790, -0.00024970776439, 0.00368937321867, -0.00466522082114, 0.00021892585197, 0.00320063483118, 0.00476416347615, 0.00157335087127;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
