#include "../include/roothaan.h"
#include "../include/system.h"
#include "../include/mp.h"

int test_numgrad_formaldehyde_mp2_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/formaldehyde.xyz", "CC-PVDZ", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-12;
    data.mp.grad.step = 1e-5, data.mp.grad.numerical = true;

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

    // transform the coulomb tensor to MO basis
    data.intsmo.J = Transform::Coulomb(data.ints.J, data.roothaan.C);

    // calculate the gradient
    libint2::initialize();
    data = MP(data).Gradient.mp2(false);
    libint2::finalize();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << -0.00102084974713, -0.00002579617836, -0.00011763843077, 0.00826570192222, 0.00020910503411, 0.00095183947947, -0.00270418815884, 0.00147472142515, -0.00873542005638, -0.00454058656549, -0.00165797541237, 0.00790121358795;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.mp.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.mp.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.mp.grad.G - G).norm() > 1e-8;
}
