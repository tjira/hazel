#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_pyrrole_hf_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/pyrrole.xyz", "CC-PVDZ", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-8;
    data.roothaan.grad.step = 0.0005, data.roothaan.grad.numerical = false;

    // initialize the guess density matrix
    data.roothaan.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

    // calculate integrals
    libint2::initialize();
    data.ints.S = Integral::Overlap(data.system);
    data.ints.T = Integral::Kinetic(data.system);
    data.ints.V = Integral::Nuclear(data.system);
    data.ints.J = Integral::Coulomb(data.system);
    data.ints.dS = Integral::dOverlap(data.system);
    data.ints.dT = Integral::dKinetic(data.system);
    data.ints.dV = Integral::dNuclear(data.system);
    data.ints.dJ = Integral::dCoulomb(data.system);
    libint2::finalize();

    // perform the SCF cycle and calculate the gradient
    data = Roothaan(data).scf(false);
    data = Roothaan(data).gradient();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << -0.00000000000000, 0.00000000000005, -0.00063134812791, -0.00000000000000, -0.01446077891416, -0.02185669490004, 0.00000000000000, 0.01446077891428, -0.02185669490005, -0.00000000000000, 0.00153569310056, 0.00863112004357, 0.00000000000000, -0.00153569310093, 0.00863112004358, 0.00000000000000, 0.00192012085373, 0.00130333837883, -0.00000000000000, -0.00192012085376, 0.00130333837885, -0.00000000000000, -0.00094334432304, 0.00010750147070, 0.00000000000000, 0.00094334432303, 0.00010750147069, -0.00000000000000, -0.00000000000000, 0.02426081814022;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
