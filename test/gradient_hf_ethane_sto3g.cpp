#include "../include/gradient.h"

int test_gradient_hf_ethane_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "STO-3G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.00000711676401, -0.00000044518034, -0.00000080732515, -0.00000713266231, -0.00000038357673, -0.00000000683923, 0.00000351306717, 0.00000117949002, 0.00000769977863, 0.00000321921503, 0.00000592083077, -0.00000504068893, 0.00000230334535, -0.00000663144425, -0.00000254525334, -0.00000314368845, -0.00000699407166, 0.00000338464372, -0.00000279589220, 0.00000637331451, 0.00000457548015, -0.00000308014877, 0.00000098063768, -0.00000725979584;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
