#include "../include/gradient.h"

int test_grad_ethylene_hf_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "CC-PVDZ", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.02004172052557, 0.00027403811039, 0.00533721082693, -0.02005653534748, -0.00028129751893, -0.00528395148981, -0.00022566318302, 0.00217754596134, 0.00026857698276, -0.00000144227393, -0.00217637892868, -0.00035720319829, 0.00000771924021, 0.00218023124955, 0.00033029735502, 0.00023420104023, -0.00217413887353, -0.00029493046853;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}