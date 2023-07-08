#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_acetone_hf_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/acetone.xyz", "CC-PVDZ", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.02044690693650, -0.00186127705101, 0.05241323652666, 0.01773381223272, 0.00114086300543, -0.01000838486964, 0.01114212979660, 0.00093648945277, -0.02855090018208, -0.00626139923978, -0.00002012259222, -0.01941044226258, 0.00048472158170, 0.00007071875268, -0.00293445929407, 0.00035805862636, -0.00300215587633, 0.00319634342479, 0.00002313172413, 0.00293539859040, 0.00328036622772, -0.00217046616460, -0.00313285902459, 0.00219766096621, -0.00249253265786, 0.00280558471181, 0.00230192741490, 0.00162945103778, 0.00012736003104, -0.00248534795527;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
