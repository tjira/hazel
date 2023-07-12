#include "../include/gradient.h"

int test_grad_ethane_hf_631gs(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "6-31G*", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00690777624291, 0.00020640064662, 0.00019134872291, 0.00690776388221, -0.00020702631294, -0.00019196277062, -0.00078870466876, 0.00005373444726, 0.00096497414054, -0.00080566210302, 0.00082569192977, -0.00047455361426, -0.00085348650540, -0.00080565975870, -0.00042216739858, 0.00080569264857, -0.00082696376574, 0.00047257549713, 0.00085310234881, 0.00080507075961, 0.00042445161234, 0.00078907064031, -0.00005124794586, -0.00096466618947;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
