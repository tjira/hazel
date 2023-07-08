#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethanol_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethanol.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00008101933320, -0.00000716672082, 0.00003789532453, -0.00004924044253, -0.00004397450753, -0.00001335930833, 0.00009993825809, 0.00003690240823, -0.00003628121678, 0.00002088059571, 0.00002877761063, -0.00000340093494, 0.00000111164037, -0.00000828048486, 0.00000644846875, -0.00000393489075, -0.00001522371730, 0.00000669204629, 0.00003341484954, 0.00002276122179, -0.00001259812765, -0.00002100578198, -0.00000526434956, 0.00001800747666, -0.00000014489839, -0.00000853146109, -0.00000340372840;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
