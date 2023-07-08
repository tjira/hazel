#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_acetone_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/acetone.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << 0.04913619885350, 0.00456108083259, -0.12597079667565, 0.04429425807482, 0.00252761932559, -0.00499763570651, -0.03317006338158, -0.00311900610361, 0.08504386975802, -0.02917015379832, -0.00106930517658, -0.03378363615094, -0.03800749177421, -0.00119820975192, -0.06475416461749, -0.01554173649729, -0.06152122401128, 0.03968364240279, -0.02217270287411, 0.05831771705750, 0.04134516096926, -0.01023891103583, -0.06132106429509, 0.04166767096417, -0.01687714733751, 0.05851536810002, 0.04350945347322, 0.07174774977119, 0.00430702402278, -0.02174356441373;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
