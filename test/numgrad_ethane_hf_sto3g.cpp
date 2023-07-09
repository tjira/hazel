#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ethane_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00000711886575, -0.00000044330760, -0.00000080351852, -0.00000713503388, -0.00000038577913, -0.00000000526404, 0.00000351412451, 0.00000117463353, 0.00000770129447, 0.00000321670609, 0.00000592505607, -0.00000503994488, 0.00000230377069, -0.00000662818177, -0.00000254666866, -0.00000314564152, -0.00000699403274, 0.00000338026742, -0.00000279370266, 0.00000637249969, 0.00000457407711, -0.00000307532895, 0.00000098136796, -0.00000726137091;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
