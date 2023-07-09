#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ethane_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "6-31G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00451031262849, 0.00013469970871, 0.00012478338027, 0.00451029420434, -0.00013531672971, -0.00012538423314, -0.00053700041225, 0.00007068556541, 0.00175471524258, -0.00056815715300, 0.00149635677872, -0.00090013626207, -0.00065605501004, -0.00151379504819, -0.00080539439751, 0.00056818309721, -0.00149767053338, 0.00089809494151, 0.00065565795082, 0.00151317915519, 0.00080773877653, 0.00053737791931, -0.00006812648865, -0.00175441218412;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
