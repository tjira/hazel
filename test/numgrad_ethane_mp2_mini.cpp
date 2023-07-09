#include "../include/roothaan.h"
#include "../include/system.h"
#include "../include/mp.h"

int test_numgrad_ethane_mp2_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.04752345007959, 0.00142115345037, 0.00131906590577, 0.04752346415179, -0.00142212086179, -0.00132015693663, -0.03282351263181, -0.00143315954759, -0.07748296462981, -0.03142619415286, -0.06574680828407, 0.04210486816116, -0.02748501670242, 0.06992443520489, 0.03792604940614, 0.03142699820196, 0.06575725833085, -0.04208799251013, 0.02748456041794, -0.06991427353850, -0.03794306401857, 0.03282314657779, 0.00141352264215, 0.07748419026830;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.mp.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.mp.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.mp.grad.G - G).norm() > 1e-8;
}
