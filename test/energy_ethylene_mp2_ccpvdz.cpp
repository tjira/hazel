#include "../include/mp.h"

int test_energy_ethylene_mp2_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "CC-PVDZ", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // calculate HF energy
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    libint2::finalize();

    // calculate the correlation energy
    double Ecorr = MP({rhfres}).rmp2(system, Tensor<4>(), false);

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << rhfres.E + Ecorr << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -78.31804852374655 << std::endl;

    // return success or failure based on the error
    return std::abs(rhfres.E + Ecorr - -78.31804852374655) > 1e-8;
}
