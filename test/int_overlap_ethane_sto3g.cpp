#include "../include/integral.h"

int test_int_overlap_ethane_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix S = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Sexp(system.shells.nbf(), system.shells.nbf()); Sexp << 1.00000000000000, 0.24836239705696, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000006653621, 0.02675582341146, -0.04481411076540, 0.00134026997636, 0.00124442068758, 0.00517610164888, 0.00517611022830, 0.00517614638720, 0.06285187814294, 0.06285188441479, 0.06285195286803, 0.24836239705696, 1.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.02675582341146, 0.29960080495282, -0.33700254184657, 0.01007884304922, 0.00935805548023, 0.09409920296933, 0.09409931495277, 0.09409978691829, 0.49302692310666, 0.49302694774478, 0.49302721665426, -0.00000000000000, -0.00000000000000, 1.00000000000000, -0.00000000000000, -0.00000000000000, 0.04481411076540, 0.33700254184657, -0.31368961202718, 0.01431880778051, 0.01329479950886, 0.10933545548439, 0.10832025201528, 0.10545692708174, -0.17049697895253, -0.14845743393106, -0.17830082526147, 0.00000000000000, 0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00134026997636, -0.01007884304922, 0.01431880778051, 0.16465482288752, -0.00039761183072, -0.00151475704136, 0.04521275762936, -0.05336175302878, -0.36776438304074, 0.39068820279384, -0.00805221195135, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00124442068758, -0.00935805548023, 0.01329479950886, -0.00039761183072, 0.16471388340605, 0.05392307669704, -0.03296657258937, -0.02992916953041, 0.23513920770535, 0.21197509394502, -0.43330641448886, 0.00000006653621, 0.02675582341146, 0.04481411076540, -0.00134026997636, -0.00124442068758, 1.00000000000000, 0.24836239705696, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.06285185733741, 0.06285187384301, 0.06285202185933, 0.00517612178011, 0.00517611133545, 0.00517612512231, 0.02675582341146, 0.29960080495282, 0.33700254184657, -0.01007884304922, -0.00935805548023, 0.24836239705696, 1.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.49302684137483, 0.49302690621496, 0.49302748767721, 0.09409946573331, 0.09409932940390, 0.09409950935755, -0.04481411076540, -0.33700254184657, -0.31368961202718, 0.01431880778051, 0.01329479950886, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.17830499730697, 0.17049276701620, 0.14845756092454, -0.10832099319172, -0.10545636484996, -0.10933528796690, 0.00134026997636, 0.01007884304922, 0.01431880778051, 0.16465482288752, -0.00039761183072, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00817197929416, 0.36770104868082, -0.39074551141203, -0.04522106813738, 0.05335409193394, 0.00153032974311, 0.00124442068758, 0.00935805548023, 0.01329479950886, -0.00039761183072, 0.16471388340605, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.43330221506214, -0.23524126781130, -0.21187005364414, 0.03295336576597, 0.02994270067996, -0.05392375910253, 0.00517610164888, 0.09409920296933, 0.10933545548439, -0.00151475704136, 0.05392307669704, 0.06285185733741, 0.49302684137483, 0.17830499730697, 0.00817197929416, 0.43330221506214, 1.00000000000000, 0.17474905680899, 0.17474915512722, 0.04765335464155, 0.04767067690901, 0.01561403053522, 0.00517611022830, 0.09409931495277, 0.10832025201528, 0.04521275762936, -0.03296657258937, 0.06285187384301, 0.49302690621496, 0.17049276701620, 0.36770104868082, -0.23524126781130, 0.17474905680899, 1.00000000000000, 0.17474933910155, 0.01561403522265, 0.04765329057675, 0.04767088576271, 0.00517614638720, 0.09409978691829, 0.10545692708174, -0.05336175302878, -0.02992916953041, 0.06285202185933, 0.49302748767721, 0.14845756092454, -0.39074551141203, -0.21187005364414, 0.17474915512722, 0.17474933910155, 1.00000000000000, 0.04767107065197, 0.01561409371910, 0.04765367777372, 0.06285187814294, 0.49302692310666, -0.17049697895253, -0.36776438304074, 0.23513920770535, 0.00517612178011, 0.09409946573331, -0.10832099319172, -0.04522106813738, 0.03295336576597, 0.04765335464155, 0.01561403522265, 0.04767107065197, 1.00000000000000, 0.17474911083511, 0.17474902478201, 0.06285188441479, 0.49302694774478, -0.14845743393106, 0.39068820279384, 0.21197509394502, 0.00517611133545, 0.09409932940390, -0.10545636484996, 0.05335409193394, 0.02994270067996, 0.04767067690901, 0.04765329057675, 0.01561409371910, 0.17474911083511, 1.00000000000000, 0.17474925650463, 0.06285195286803, 0.49302721665426, -0.17830082526147, -0.00805221195135, -0.43330641448886, 0.00517612512231, 0.09409950935755, -0.10933528796690, 0.00153032974311, -0.05392375910253, 0.01561403053522, 0.04767088576271, 0.04765367777372, 0.17474902478201, 0.17474925650463, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP: " << S << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP NORM: " << S.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP: " << Sexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP NORM: " << Sexp.norm() << std::endl;

    // return success or failure based on the error
    return (S - Sexp).norm() > 1e-8;
}
