#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_ethanol_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethanol.xyz", "3-21G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00663468744348, 0.01936999952811, -0.02486727886012, -0.00508391566304, -0.00131434343491, 0.00129582120629, 0.01341645852321, -0.02074531312150, 0.00611512185557, -0.00131056157107, 0.00084884994193, -0.00043198047432, 0.00028583743081, -0.00034300898919, 0.00030912705316, -0.00028197527497, 0.00416319024438, 0.00653880762014, -0.01128562703945, -0.00767931778694, 0.01961510304816, -0.00079291079194, 0.00683696511835, -0.00678819349082, -0.00158199305757, -0.00113702149955, -0.00178652795835;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
