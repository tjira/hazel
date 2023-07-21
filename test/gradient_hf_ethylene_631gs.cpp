#include "../include/gradient.h"

int test_gradient_hf_ethylene_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "6-31G*", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << 0.01986478570431, 0.00027130545157, 0.00529303847448, -0.01988097524543, -0.00027921361421, -0.00523476656546, -0.00264945793781, -0.00315046557715, -0.00120737025915, -0.00297468725403, 0.00307703557060, -0.00031352868692, 0.00298152539331, -0.00307284451785, 0.00028408465456, 0.00265880933928, 0.00315418268703, 0.00117854238261;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
