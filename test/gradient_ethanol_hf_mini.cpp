#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethanol_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethanol.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.03675660549682, -0.01105738589604, 0.05083061677561, -0.03805267932631, -0.02346278426006, -0.00320297959466, 0.00612142428990, 0.04124016549786, 0.00556810231282, 0.05869595172052, -0.04525333397854, -0.00603146864807, 0.00499145931716, 0.04866778807369, -0.05469945524447, 0.00054889163591, -0.03954169694707, -0.06245615669754, -0.00354225768497, 0.04551480546319, -0.04852437572780, 0.00285577701903, -0.05296899202988, 0.05432498031601, 0.00513803852656, 0.03686143407682, 0.06419073650775;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
