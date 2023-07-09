#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ammonia_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ammonia.xyz", "MINI", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-12;
    data.roothaan.grad.step = 1e-5, data.roothaan.grad.numerical = true;

    // initialize the guess density matrix
    data.roothaan.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // calculate integrals
    libint2::initialize();
    data.ints.S = Integral::Overlap(data.system);
    data.ints.T = Integral::Kinetic(data.system);
    data.ints.V = Integral::Nuclear(data.system);
    data.ints.J = Integral::Coulomb(data.system);
    libint2::finalize();

    // perform the SCF cycle
    data = Roothaan(data).scf(false);

    // calculate the gradient
    libint2::initialize();
    data = Roothaan(data).gradient(false);
    libint2::finalize();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << 0.03261452450307, 0.04060336749805, -0.02632376843529, -0.01714937259489, -0.04094025401093, -0.04127738417379, -0.04859921956588, 0.02684709230059, 0.02431857180476, 0.03313406822171, -0.02651020183967, 0.04328258362434;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
