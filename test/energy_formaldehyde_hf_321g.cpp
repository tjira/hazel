#include "../include/hf.h"

int test_energy_formaldehyde_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/formaldehyde.xyz", "3-21G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // perform the SCF cycle
    libint2::initialize();
    data = HF(data).rscf(false);
    libint2::finalize();

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << data.hf.E << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -113.22119680543543 << std::endl;

    // return success or failure based on the error
    return std::abs(data.hf.E - -113.22119680543543) > 1e-8;
}
