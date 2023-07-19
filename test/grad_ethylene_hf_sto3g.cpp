#include "../include/gradient.h"

int test_grad_ethylene_hf_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "STO-3G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.00005808196918, -0.00000404350559, 0.00004935715373, -0.00007580128277, -0.00000459825379, 0.00001408672869, 0.00001981100387, 0.00003584827035, -0.00000651535329, 0.00002109561221, -0.00003009510185, -0.00001656575247, -0.00001354558672, 0.00003466338865, -0.00001547928935, -0.00000964171441, -0.00003177479757, -0.00002488348761;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
