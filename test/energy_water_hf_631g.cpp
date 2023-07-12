#include "../include/hf.h"

int test_energy_water_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/water.xyz", "6-31G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // perform the SCF cycle
    libint2::initialize();
    data = HF(data).scf(false);
    libint2::finalize();

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << data.hf.E << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -75.97974730915679 << std::endl;

    // return success or failure based on the error
    return std::abs(data.hf.E - -75.97974730915679) > 1e-8;
}
