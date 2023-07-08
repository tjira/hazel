#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_acetone_hf_sto3g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/acetone.xyz", "STO-3G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.00003085070516, 0.00002819544257, -0.00008377465004, -0.00000614562882, 0.00002681903907, -0.00012649034361, -0.00009018403136, -0.00008600136463, 0.00024257609128, 0.00008817389623, 0.00003178355568, -0.00009021689993, 0.00002335522425, 0.00000199081651, 0.00003888872053, -0.00000478951768, -0.00000118662509, -0.00000437221014, -0.00001460712087, -0.00000101895269, 0.00000019177397, 0.00000527496892, 0.00000009475532, 0.00000104778304, 0.00001156912519, 0.00000170772599, 0.00000879365357, -0.00004349762167, -0.00000238439266, 0.00001335608128;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
