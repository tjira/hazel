#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ethane_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "3-21G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00303643444012, -0.00009100401824, -0.00008477439971, -0.00303644271219, 0.00009039038126, 0.00008418445092, 0.00013478882143, 0.00004068729006, 0.00141416890400, 0.00010936161516, 0.00120227125705, -0.00074967112147, 0.00003764204136, -0.00125085686718, -0.00067205995587, -0.00010931386277, -0.00120355756348, 0.00074766589720, -0.00003804361262, 0.00125024924625, 0.00067435658251, -0.00013439965827, -0.00003816731762, -0.00141386396553;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
