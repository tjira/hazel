#include "../include/mp.h"

int test_energy_ethylene_mp2_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "3-21G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // calculate HF energy and MP2 correlation
    libint2::initialize();
    data = MP(HF(data).scf(false)).mp2(false);
    libint2::finalize();

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << data.hf.E + data.mp.Ecorr << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -77.78309713802911 << std::endl;

    // return success or failure based on the error
    return std::abs(data.hf.E + data.mp.Ecorr - -77.78309713802911) > 1e-8;
}
