#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethane_hf_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "CC-PVDZ", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.01298130182916, 0.00038805275048, 0.00035999774064, 0.01298129236340, -0.00038869508382, -0.00036064351804, -0.00267587955401, -0.00004113722206, -0.00389560793700, -0.00260526492222, -0.00329975616725, 0.00215941898432, -0.00240606271088, 0.00357137825883, 0.00194995582081, 0.00260533310249, 0.00329910978844, -0.00216039419700, 0.00240567839905, -0.00357139039212, -0.00194869562320, 0.00267620515134, 0.00004243806749, 0.00389596872946;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
