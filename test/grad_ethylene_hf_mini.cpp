#include "../include/gradient.h"

int test_grad_ethylene_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.02623301087941, 0.00036299197068, 0.00694641059751, -0.02623234386978, -0.00036295331880, -0.00695041385935, 0.03626700289653, 0.06925390458595, 0.02041880184122, 0.04340009910297, -0.06814744824546, 0.00068585669171, -0.04339977047471, 0.06814745648944, -0.00068368080703, -0.03626799853802, -0.06925395148226, -0.02041697446727;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
