#include "../include/integral.h"

int test_integral_overlap_water_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.23368985719701, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.16727976258450, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.03353153616877, -0.00000000000000, 0.00000000000000, 0.03353153616877, -0.00000000000000, 0.03353153616877, 0.03064514778093, 0.06662563164355, 0.03064494066984, 0.06662551280329, 0.23368985719701, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.76364080963424, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.54706575927831, 0.00000000000000, 0.00000000000000, 0.54706575927831, 0.00000000000000, 0.54706575927831, 0.22740089882644, 0.36894290363964, 0.22739985277649, 0.36894233001181, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.18911105089794, 0.10411698434865, 0.11627108451705, 0.06401443574809, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.13629852194184, -0.07504051724308, 0.24209657702555, 0.13328916505085, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.13389599191494, 0.07371778025853, -0.01175548476966, -0.00647212269156, 0.16727976258450, 0.76364080963424, -0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.69901474023950, 0.40941140833259, 0.66883637172093, 0.40941050365403, 0.66883561035136, -0.00000000000000, -0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.40622262390300, 0.34187322277300, 0.24975843830849, 0.21019463272268, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.29277793632059, -0.24639921746615, 0.52004041459565, 0.43766170499459, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.28761714827091, 0.24205594572604, -0.02525160515897, -0.02125154172166, 0.03353153616877, 0.54706575927831, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.33333333333333, 0.34180652199309, 0.42469048752734, 0.22519840942041, 0.40207743218542, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.13512090694470, -0.02620263352568, 0.14756292957896, 0.02861547993432, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.13273913470262, 0.02574076047722, -0.00716521395116, -0.00138948201035, 0.03353153616877, 0.54706575927831, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.33333333333333, 0.25171544338367, 0.40722003406375, 0.46158037963137, 0.44791674660338, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, -0.09566943749665, -0.01855220829267, -0.01491921898241, -0.00289314269271, 0.03353153616877, 0.54706575927831, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 1.00000000000000, 0.24831245984563, 0.40656012781713, 0.15505318302448, 0.38847483345586, 0.03064514778093, 0.22740089882644, 0.18911105089794, -0.13629852194184, 0.13389599191494, 0.40941140833259, 0.40622262390300, -0.29277793632059, 0.28761714827091, 0.34180652199309, -0.13512090694470, 0.13273913470262, 0.25171544338367, -0.09566943749665, 0.24831245984563, 1.00000000000000, 0.65829196968307, 0.05152895979194, 0.22404505366272, 0.06662563164355, 0.36894290363964, 0.10411698434865, -0.07504051724308, 0.07371778025853, 0.66883637172093, 0.34187322277300, -0.24639921746615, 0.24205594572604, 0.42469048752734, -0.02620263352568, 0.02574076047722, 0.40722003406375, -0.01855220829267, 0.40656012781713, 0.65829196968307, 1.00000000000000, 0.22404505366272, 0.51583577082327, 0.03064494066984, 0.22739985277649, 0.11627108451705, 0.24209657702555, -0.01175548476966, 0.40941050365403, 0.24975843830849, 0.52004041459565, -0.02525160515897, 0.22519840942041, 0.14756292957896, -0.00716521395116, 0.46158037963137, -0.01491921898241, 0.15505318302448, 0.05152895979194, 0.22404505366272, 1.00000000000000, 0.65829196968307, 0.06662551280329, 0.36894233001181, 0.06401443574809, 0.13328916505085, -0.00647212269156, 0.66883561035136, 0.21019463272268, 0.43766170499459, -0.02125154172166, 0.40207743218542, 0.02861547993432, -0.00138948201035, 0.44791674660338, -0.00289314269271, 0.38847483345586, 0.22404505366272, 0.51583577082327, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}