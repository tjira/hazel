#include "../include/gradient.h"

int test_grad_formaldehyde_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/formaldehyde.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00000840308643, 0.00000023378289, 0.00000085522386, -0.00000710463653, -0.00000021162035, -0.00000065428938, -0.00000037260803, 0.00000036430138, -0.00000207964827, -0.00000092584186, -0.00000038646387, 0.00000187871382;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
