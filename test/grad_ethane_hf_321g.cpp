#include "../include/gradient.h"

int test_grad_ethane_hf_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "3-21G", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // calculate HF energy and gradient
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    Matrix G = Gradient({}).get(system, rhfres, false);
    libint2::finalize();

    // create the expectation gradient
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.00303642532834, -0.00009100993693, -0.00008478441460, -0.00303644033618, 0.00009039331134, 0.00008418305304, 0.00013477102236, 0.00004069524851, 0.00141418036805, 0.00010934789090, 0.00120227811305, -0.00074967357460, 0.00003763016835, -0.00125085217992, -0.00067205713948, -0.00010930917666, -0.00120356950233, 0.00074767221515, -0.00003803385953, 0.00125024624512, 0.00067435733951, -0.00013439103683, -0.00003818129885, -0.00141387784707;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
