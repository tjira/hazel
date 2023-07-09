#include "../include/roothaan.h"
#include "../include/system.h"

int test_numgrad_ethane_hf_mini(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "MINI", 0, 1);

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
    Matrix G(data.system.atoms.size(), 3); G << -0.04439093890753, 0.00132745395255, 0.00123206310672, 0.04439094943561, -0.00132842028039, -0.00123313659543, -0.02866209812825, -0.00125347234440, -0.06772838978526, -0.02744072554025, -0.05746998241515, 0.03680283187940, -0.02399583372265, 0.06111991020004, 0.03315036870594, 0.02744144220207, 0.05747900761637, -0.03678824897710, 0.02399536334284, -0.06111102524797, -0.03316502004079, 0.02866183718214, 0.00123653566301, 0.06772952719448;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << data.roothaan.grad.G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << data.roothaan.grad.G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << G.norm() << std::endl;

    // return success or failure based on the error
    return (data.roothaan.grad.G - G).norm() > 1e-8;
}
