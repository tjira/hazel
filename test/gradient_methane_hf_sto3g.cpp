#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_methane_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/methane.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00000144027520, 0.00000124042828, 0.00000196539748, -0.00000013270742, -0.00000059119857, -0.00000004036382, -0.00000073796858, -0.00000070037504, -0.00000073527914, -0.00000036693637, 0.00000016687391, -0.00000042516860, -0.00000020266282, -0.00000011572858, -0.00000076458592;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
