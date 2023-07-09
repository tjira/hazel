#include "../include/roothaan.h"
#include "../include/system.h"
#include "../include/mp.h"

int test_numgrad_ethylene_mp2_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethylene.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.04047531822872, 0.00055785820032, 0.01074533358422, -0.04048808105191, -0.00056407455814, -0.01070021215775, 0.00623838027640, 0.01213576407541, 0.00353485929494, 0.00748658947110, -0.01194134924731, 0.00007654146556, -0.00748117267386, 0.01194466650778, -0.00009933580924, -0.00623103585362, -0.01213287090085, -0.00355718527616;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.mp.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.mp.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.mp.grad.G - G).norm() > 1e-8;
}
