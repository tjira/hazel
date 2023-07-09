#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ethane_hf_ccpvdz(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "CC-PVDZ", 0, 1);

    // set some options
    data.roothaan.diis = {3, 5}, data.roothaan.maxiter = 1000, data.roothaan.thresh = 1e-12;
    data.roothaan.grad.step = 1e-5, data.roothaan.grad.numerical = true;

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

    // calculate the gradient
    libint2::initialize();
    data = Roothaan(data).gradient(false);
    libint2::finalize();

    // create the expectation gradient
    Matrix G(data.system.atoms.size(), 3); G << -0.01298150192435, 0.00038806297204, 0.00036000224040, 0.01298150154835, -0.00038869315315, -0.00036064558162, -0.00267588661929, -0.00004113172567, -0.00389560028183, -0.00260526798242, -0.00329975764825, 0.00215941860560, -0.00240606795847, 0.00357137322764, 0.00194996460247, 0.00260533829499, 0.00329910753897, -0.00216039846955, 0.00240567841931, -0.00357139729184, -0.00194869559218, 0.00267619870182, 0.00004243570425, 0.00389596989283;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
