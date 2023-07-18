#include "../include/gradient.h"

int test_grad_methane_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/methane.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00000144027520, 0.00000124042828, 0.00000196539748, -0.00000013270742, -0.00000059119857, -0.00000004036382, -0.00000073796858, -0.00000070037504, -0.00000073527914, -0.00000036693637, 0.00000016687391, -0.00000042516860, -0.00000020266282, -0.00000011572858, -0.00000076458592;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
