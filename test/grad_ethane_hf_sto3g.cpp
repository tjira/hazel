#include "../include/roothaan.h"
#include "../include/system.h"

int test_grad_ethane_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "STO-3G", 0, 1);

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

    // perform the SCF cycle
    data = Roothaan(data).scf(false);

    // calculate the gradient
    data = Roothaan(data).gradient(false);

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << 0.00000711676384, -0.00000044518035, -0.00000080732515, -0.00000713266195, -0.00000038357673, -0.00000000683924, 0.00000351306717, 0.00000117949002, 0.00000769977863, 0.00000321921504, 0.00000592083077, -0.00000504068893, 0.00000230334536, -0.00000663144425, -0.00000254525334, -0.00000314368845, -0.00000699407166, 0.00000338464372, -0.00000279589220, 0.00000637331451, 0.00000457548015, -0.00000308014877, 0.00000098063768, -0.00000725979585;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
