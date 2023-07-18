#include "../include/mp.h"

int test_energy_water_mp2_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/water.xyz", "STO-3G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // calculate HF energy and MP2 correlation
    libint2::initialize();
    data = MP(HF(data).rscf(false)).mp2(false);
    libint2::finalize();

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << data.hf.E + data.mp.Ecorr << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -75.00485494802203 << std::endl;

    // return success or failure based on the error
    return std::abs(data.hf.E + data.mp.Ecorr - -75.00485494802203) > 1e-8;
}