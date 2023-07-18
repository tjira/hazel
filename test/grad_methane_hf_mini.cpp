#include "../include/gradient.h"

int test_grad_methane_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/methane.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00000065169929, 0.00000050858133, 0.00000064780829, -0.03746567190673, -0.01141011005943, -0.06630736208987, -0.01453502144294, -0.05838127098654, 0.04807158028413, -0.02370089945815, 0.06598735460504, 0.03184997217409, 0.07570094110853, 0.00380351785960, -0.01361483817664;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
