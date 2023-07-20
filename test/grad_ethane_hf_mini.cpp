#include "../include/gradient.h"

int test_grad_ethane_hf_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "MINI", 0, 1);

    // set some options
    HF::OptionsRestricted rhfopt = {{3, 5}, 1e-8, 1000, false};

    // initialize the guess density matrix
    Matrix D(system.shells.nbf(), system.shells.nbf());

    // calculate HF energy and gradient
    libint2::initialize();
    HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
    Matrix G = Gradient({}).get(system, rhfres, false);
    libint2::finalize();

    // create the expectation gradient
    Matrix Gexp(system.atoms.size(), 3); Gexp << -0.04439093638970, 0.00132745202383, 0.00123206536890, 0.04439095001374, -0.00132842262872, -0.00123313349451, -0.02866209762443, -0.00125347131433, -0.06772839404478, -0.02744072771361, -0.05746998438062, 0.03680283036100, -0.02399583308136, 0.06111990759537, 0.03315036968817, 0.02744144367081, 0.05747901062801, -0.03678824650985, 0.02399535946431, -0.06111102499627, -0.03316501962251, 0.02866184166025, 0.00123653307273, 0.06772952825360;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
