#include "../include/roothaan.h"
#include "../include/system.h"
#include "../include/mp.h"

int test_numgrad_ammonia_mp2_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ammonia.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.04810468542539, 0.05988792700425, -0.03882584932320, -0.02361333402326, -0.05304595931177, -0.04747895863080, -0.06157848599200, 0.02878456462326, 0.03170605048590, 0.03708713536025, -0.03562652820889, 0.05459876017503;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.mp.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.mp.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.mp.grad.G - G).norm() > 1e-8;
}
