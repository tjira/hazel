#include "../include/gradient.h"

int test_grad_ethylene_hf_631gs(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "6-31G*", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.01986478570408, 0.00027130545134, 0.00529303847403, -0.01988097524595, -0.00027921361475, -0.00523476655702, -0.00264945793780, -0.00315046557715, -0.00120737025913, -0.00297468725402, 0.00307703557060, -0.00031352868691, 0.00298152539332, -0.00307284451785, 0.00028408465460, 0.00265880933929, 0.00315418268702, 0.00117854238265;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
