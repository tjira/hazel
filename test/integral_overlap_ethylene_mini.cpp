#include "../include/integral.h"

int test_integral_overlap_ethylene_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.19254595363088, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000013553934, 0.04732025011986, -0.08592728331793, -0.00119212438465, -0.02276058834275, 0.06909899354342, 0.06909882734529, 0.01022010646639, 0.01021986792309, 0.19254595363088, 1.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.04732025011986, 0.46908373490887, -0.46810314634570, -0.00649429556877, -0.12399208498770, 0.58084214420611, 0.58084144884479, 0.17567697580744, 0.17567434843876, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.08592728331793, 0.46810314634570, -0.20173964437270, -0.00766496101476, -0.14634296938067, -0.23831565322149, -0.28225770863001, 0.21033261399645, 0.19973813282325, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00119212438465, 0.00649429556877, -0.00766496101476, 0.35063770282845, -0.00203030999682, -0.42671906764612, 0.41951708669072, -0.09917636587982, 0.10486883255056, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.02276058834275, 0.12399208498770, -0.14634296938067, -0.00203030999682, 0.31198042928583, -0.12980426757552, -0.00827455422550, 0.03963554423089, 0.06893915638410, 0.00000013553934, 0.04732025011986, 0.08592728331793, 0.00119212438465, 0.02276058834275, 1.00000000000000, 0.19254595363088, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.01022019300581, 0.01021985419714, 0.06909880748777, 0.06909898684267, 0.04732025011986, 0.46908373490887, 0.46810314634570, 0.00649429556877, 0.12399208498770, 0.19254595363088, 1.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.17567792896732, 0.17567419725744, 0.58084136576222, 0.58084211617066, -0.08592728331793, -0.46810314634570, -0.20173964437270, -0.00766496101476, -0.14634296938067, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.19972909688805, -0.21031962526478, 0.28229580661886, 0.23837234563833, -0.00119212438465, -0.00649429556877, -0.00766496101476, 0.35063770282845, -0.00203030999682, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.10486620152594, 0.09917990990906, -0.41949489823803, 0.42674005011664, -0.02276058834275, -0.12399208498770, -0.14634296938067, -0.00203030999682, 0.31198042928583, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.06898252329091, -0.03967776311996, 0.00809649535394, 0.12963104501800, 0.06909899354342, 0.58084214420611, -0.23831565322149, -0.42671906764612, -0.12980426757552, 0.01022019300581, 0.17567792896732, -0.19972909688805, -0.10486620152594, -0.06898252329091, 1.00000000000000, 0.25676436340716, 0.10753177232239, 0.03881851850550, 0.06909882734529, 0.58084144884479, -0.28225770863001, 0.41951708669072, -0.00827455422550, 0.01021985419714, 0.17567419725744, -0.21031962526478, 0.09917990990906, -0.03967776311996, 0.25676436340716, 1.00000000000000, 0.03881823789753, 0.10752499477268, 0.01022010646639, 0.17567697580744, 0.21033261399645, -0.09917636587982, 0.03963554423089, 0.06909880748777, 0.58084136576222, 0.28229580661886, -0.41949489823803, 0.00809649535394, 0.10753177232239, 0.03881823789753, 1.00000000000000, 0.25676528425895, 0.01021986792309, 0.17567434843876, 0.19973813282325, 0.10486883255056, 0.06893915638410, 0.06909898684267, 0.58084211617066, 0.23837234563833, 0.42674005011664, 0.12963104501800, 0.03881851850550, 0.10752499477268, 0.25676528425895, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}