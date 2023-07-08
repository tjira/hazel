#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethanol_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethanol.xyz", "6-31G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.01520079772765, 0.02465873724342, -0.04125403550683, -0.01276455469639, -0.00385753284856, 0.00141419719193, 0.01867701819618, -0.01808713824780, 0.01074191613773, -0.00102886641049, 0.00114483189879, -0.00044435428224, 0.00101688391697, -0.00020745368533, 0.00049743685505, -0.00289584370264, 0.00506412854229, 0.00680486861278, -0.01568293209501, -0.01481220887918, 0.03161061813700, -0.00201517513701, 0.00713444149969, -0.00735385199550, -0.00050732779981, -0.00103780552373, -0.00201679514998;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
