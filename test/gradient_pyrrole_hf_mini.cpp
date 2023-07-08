#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_pyrrole_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/pyrrole.xyz", "MINI", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-8;
    data.roothaan.grad.step = 0.0005, data.roothaan.grad.numerical = false;

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00000000000000, 0.00000000000003, 0.00882293918841, -0.00000000000000, -0.01512811129448, -0.02552453330991, 0.00000000000000, 0.01512811129447, -0.02552453330991, -0.00000000000000, 0.01479580914584, 0.01906105773319, 0.00000000000000, -0.01479580914568, 0.01906105773319, -0.00000000000000, -0.07457130236893, -0.02937036355280, 0.00000000000000, 0.07457130236892, -0.02937036355279, 0.00000000000000, -0.04584168284381, 0.06285292448344, 0.00000000000000, 0.04584168284380, 0.06285292448343, -0.00000000000000, 0.00000000000000, -0.06286110989675;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
