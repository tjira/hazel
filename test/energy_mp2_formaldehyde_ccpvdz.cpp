#include "../include/mp.h"

int test_energy_mp2_formaldehyde_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "CC-PVDZ", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // perform the SCF cycle
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    libint2::finalize();

    // calculate the correlation energy
    double Ecorr = MP({rhfres}).rmp2(system, Tensor<4>(), false);

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << rhfres.E + Ecorr << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -114.19754586038211 << std::endl;

    // return success or failure based on the error
    return std::abs(rhfres.E + Ecorr - -114.19754586038211) > 1e-8;
}
