#include "../include/roothaan.h"
#include "../include/system.h"
#include "../include/mp.h"

int test_numgrad_ethane_mp2_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "CC-PVDZ", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.01636856479155, 0.00048937855815, 0.00045406434874, 0.01636856829455, -0.00049002016182, -0.00045473931680, -0.00481006008488, -0.00015221978821, -0.00951224124478, -0.00463825832954, -0.00806715049031, 0.00520122731901, -0.00415365853222, 0.00862676944536, 0.00468907075311, 0.00463837552205, 0.00806733921594, -0.00520086489319, 0.00415328098109, -0.00862603485290, -0.00468917570203, 0.00481030885349, 0.00015193727626, 0.00951267372760;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.mp.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.mp.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.mp.grad.G - G).norm() > 1e-8;
}
