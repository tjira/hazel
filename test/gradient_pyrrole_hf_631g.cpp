#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_pyrrole_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/pyrrole.xyz", "6-31G", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-8;

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00000000000000, -0.00000000000003, -0.01437654409476, -0.00000000000000, 0.00195202801868, -0.02287094398967, 0.00000000000000, -0.00195202801935, -0.02287094398967, -0.00000000000000, -0.00259971658351, 0.01592536731174, 0.00000000000000, 0.00259971658331, 0.01592536731171, 0.00000000000000, 0.00908874894853, 0.00503583711903, -0.00000000000000, -0.00908874894854, 0.00503583711903, -0.00000000000000, 0.00412392413759, -0.00588007990893, -0.00000000000000, -0.00412392413760, -0.00588007990893, -0.00000000000000, -0.00000000000000, 0.02995618302939;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
