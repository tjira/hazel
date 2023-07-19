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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.02004172052557, 0.00027403811038, 0.00533721082693, -0.02005653534748, -0.00028129751893, -0.00528395148981, -0.00022566318302, 0.00217754596134, 0.00026857698276, -0.00000144227393, -0.00217637892868, -0.00035720319829, 0.00000771924020, 0.00218023124956, 0.00033029735502, 0.00023420104023, -0.00217413887353, -0.00029493046853;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
