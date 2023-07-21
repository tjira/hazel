#include "../include/integral.h"

int test_integral_nuclear_water_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -61.68814245885452, -7.43627733601170, -0.01650498869282, -0.00571807049929, -0.00660132732899, -1.61453885098809, -1.61453261795134, -7.43627733601170, -10.11005794527141, -0.19630389578691, -0.06800856538633, -0.07851359260400, -3.65607297392040, -3.65606289954745, -0.01650498869282, -0.19630389578691, -10.03853985047514, -0.00380849453494, -0.03843912310759, -1.98305336034885, -1.24398827561960, -0.00571807049929, -0.06800856538633, -0.00380849453494, -10.08332042553974, 0.03385200617768, 1.36068628763617, -2.47868101188753, -0.00660132732899, -0.07851359260400, -0.03843912310759, 0.03385200617768, -9.98844868638444, -1.38426493580778, 0.09357959412604, -1.61453885098809, -3.65607297392040, -1.98305336034885, 1.36068628763617, -1.38426493580778, -5.71597436990229, -1.56781922346295, -1.61453261795134, -3.65606289954745, -1.24398827561960, -2.47868101188753, 0.09357959412604, -1.56781922346295, -5.71596855040002;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}
