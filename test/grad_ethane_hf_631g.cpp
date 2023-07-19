#include "../include/gradient.h"

int test_grad_ethane_hf_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "6-31G", 0, 1);

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
    Matrix Gexp(system.atoms.size(), 3); Gexp << -0.00451030213895, 0.00013469585995, 0.00012477760383, 0.00451028723375, -0.00013530637920, -0.00012537614977, -0.00053699169423, 0.00007068042874, 0.00175470492798, -0.00056815099756, 0.00149634628825, -0.00090014038261, -0.00065604552206, -0.00151380980068, -0.00080538750751, 0.00056818311958, -0.00149766527780, 0.00089809086235, 0.00065565058737, 0.00151318284218, 0.00080773193611, 0.00053736941190, -0.00006812396143, -0.00175440129040;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT: " << G << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED GRADIENT NORM: " << G.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT: " << Gexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED GRADIENT NORM: " << Gexp.norm() << std::endl;

    // return success or failure based on the error
    return (G - Gexp).norm() > 1e-8;
}
