#include "../include/hf.h"

int test_energy_ethane_hf_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "CC-PVDZ", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // perform the SCF cycle
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    libint2::finalize();

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << rhfres.E << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -79.23464790630865 << std::endl;

    // return success or failure based on the error
    return std::abs(rhfres.E - -79.23464790630865) > 1e-8;
}
