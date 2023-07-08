#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_acetone_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/acetone.xyz", "6-31G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00277102089144, -0.00023003838990, 0.00709907823059, 0.01505900880173, 0.00097084917569, -0.00832304500857, -0.00688875038090, -0.00072735702683, 0.01767317250371, -0.00543615954223, -0.00002209536666, -0.01635380233949, 0.00313147535932, 0.00014098481869, 0.00250859365541, 0.00082432354690, 0.00171233994517, -0.00028968482439, 0.00099895090229, -0.00160087853367, -0.00033236937307, -0.00055462600173, 0.00164293431740, -0.00082845919903, -0.00036449361424, -0.00166900722361, -0.00086849197674, -0.00399870817981, -0.00021773171633, -0.00028499166663;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
