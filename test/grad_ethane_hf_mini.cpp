#include "../include/gradient.h"

int test_grad_ethane_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.04439093638979, 0.00132745202381, 0.00123206536889, 0.04439095001330, -0.00132842262872, -0.00123313349452, -0.02866209762444, -0.00125347131433, -0.06772839404478, -0.02744072771362, -0.05746998438062, 0.03680283036100, -0.02399583308137, 0.06111990759537, 0.03315036968817, 0.02744144367081, 0.05747901062801, -0.03678824650985, 0.02399535946430, -0.06111102499627, -0.03316501962251, 0.02866184166024, 0.00123653307273, 0.06772952825360;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.hf.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.hf.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.hf.grad.G - G).norm() > 1e-8;
}