#include "../include/roothaan.h"
#include "../include/system.h"

int test_energy_ethane_hf_631g(int, char**) {
    // initialize the system
    Data data; data.system = System("../example/molecule/ethane.xyz", "6-31G", 0, 1);

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
    libint2::finalize();

    // perform the SCF cycle
    data = Roothaan(data).scf(false);

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED ENERGY: " << data.roothaan.E << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED ENERGY: " << -79.19748091997050 << std::endl;

    // return success or failure based on the error
    return std::abs(data.roothaan.E - -79.19748091997050) > 1e-8;
}
