#include "../include/gradient.h"

int test_grad_ethylene_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "3-21G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;
    data.hf.grad.step = 0.0005, data.hf.grad.numerical = false;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // perform the SCF cycle and calculate gradient
    libint2::initialize();
    data = Gradient<HF>(HF(data).rscf(false)).get(false);
    libint2::finalize();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << 0.01819013524636, 0.00024787079493, 0.00485083129517, -0.01820703117586, -0.00025627949988, -0.00479000652243, -0.00319069635535, -0.00476019863803, -0.00160339373163, -0.00368241576347, 0.00466966468791, -0.00024970776440, 0.00368937321867, -0.00466522082114, 0.00021892585196, 0.00320063483119, 0.00476416347616, 0.00157335087126;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
