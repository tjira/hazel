#include "../include/roothaan.h"
#include "../include/system.h"

int test_gradient_acetone_hf_321g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/acetone.xyz", "3-21G", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.00622048417083, -0.00054699001396, 0.01594211022879, 0.00893142049430, 0.00054749295367, -0.00263059701724, 0.00155495804854, 0.00005451191112, -0.00397383051820, -0.00478538542223, -0.00010806565138, -0.00800512293882, 0.00391588545063, 0.00019500930190, 0.00178418518095, 0.00126527469735, 0.00121506626144, 0.00002439502212, 0.00138379271958, -0.00106671469212, -0.00000414181492, -0.00104702137342, 0.00109858123695, -0.00088108279455, -0.00091351151465, -0.00118195472527, -0.00090533036726, -0.00408492892777, -0.00020693658238, -0.00135058497822;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
