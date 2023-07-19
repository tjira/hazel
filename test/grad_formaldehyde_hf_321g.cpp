#include "../include/gradient.h"

int test_grad_formaldehyde_hf_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "3-21G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << -0.01958255288120, -0.00049539947789, -0.00225484819105, 0.00483343038054, 0.00012225003361, 0.00055668084847, 0.00630259733685, -0.00164229325287, 0.01056149889293, 0.00844652516382, 0.00201544269717, -0.00886333155048;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
