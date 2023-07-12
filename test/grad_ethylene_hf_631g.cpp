#include "../include/gradient.h"

int test_grad_ethylene_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "6-31G", 0, 1);

    // set some options
    data.hf.diis = {3, 5}, data.hf.maxiter = 1000, data.hf.thresh = 1e-8;
    data.hf.grad.step = 0.0005, data.hf.grad.numerical = false;

    // initialize the guess density matrix
    data.hf.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // perform the SCF cycle and calculate gradient
    libint2::initialize();
    data = Gradient<HF>(HF(data).scf(false)).get(false);
    libint2::finalize();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << 0.02798173157832, 0.00038367453843, 0.00744479691621, -0.02799879129733, -0.00039214202258, -0.00738331496417, -0.00277881355682, -0.00471290456436, -0.00148793227535, -0.00326616819105, 0.00463389844813, -0.00014616106593, 0.00327326466556, -0.00462942466656, 0.00011506577323, 0.00278877680022, 0.00471689826615, 0.00145754562424;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
