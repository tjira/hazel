#include "../include/integral.h"

int test_int_overlap_water_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix S = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Sexp(system.shells.nbf(), system.shells.nbf()); Sexp << 1.00000000000000, 0.23670392057273, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.04999795173196, 0.04999776172299, 0.23670392057273, 1.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.45391040352827, 0.45390938059558, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.26875188293754, 0.16523695845837, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.19369827537600, 0.34405202468906, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.19028396158799, -0.01670613598051, 0.04999795173196, 0.45391040352827, 0.26875188293754, -0.19369827537600, 0.19028396158799, 1.00000000000000, 0.25089795240148, 0.04999776172299, 0.45390938059558, 0.16523695845837, 0.34405202468906, -0.01670613598051, 0.25089795240148, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP: " << S << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP NORM: " << S.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP: " << Sexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP NORM: " << Sexp.norm() << std::endl;

    // return success or failure based on the error
    return (S - Sexp).norm() > 1e-8;
}
