#include "../include/gradient.h"

int test_grad_ethylene_hf_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "CC-PVDZ", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.02004172052247, 0.00027403811033, 0.00533721082410, -0.02005653534573, -0.00028129751898, -0.00528395149497, -0.00022566318304, 0.00217754596134, 0.00026857698275, -0.00000144227394, -0.00217637892869, -0.00035720319830, 0.00000771924021, 0.00218023124955, 0.00033029735497, 0.00023420104024, -0.00217413887354, -0.00029493046858;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
