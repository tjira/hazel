#include "../include/gradient.h"

int test_grad_ethane_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "6-31G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00451030213895, 0.00013469585995, 0.00012477760383, 0.00451028723376, -0.00013530637920, -0.00012537614977, -0.00053699169423, 0.00007068042874, 0.00175470492798, -0.00056815099756, 0.00149634628825, -0.00090014038261, -0.00065604552206, -0.00151380980068, -0.00080538750751, 0.00056818311958, -0.00149766527780, 0.00089809086235, 0.00065565058737, 0.00151318284219, 0.00080773193612, 0.00053736941189, -0.00006812396143, -0.00175440129040;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}
