#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethane_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "3-21G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00303642532836, -0.00009100993692, -0.00008478441462, -0.00303644033621, 0.00009039331132, 0.00008418305305, 0.00013477102236, 0.00004069524851, 0.00141418036806, 0.00010934789090, 0.00120227811305, -0.00074967357460, 0.00003763016835, -0.00125085217992, -0.00067205713948, -0.00010930917666, -0.00120356950233, 0.00074767221515, -0.00003803385953, 0.00125024624512, 0.00067435733951, -0.00013439103682, -0.00003818129885, -0.00141387784708;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
