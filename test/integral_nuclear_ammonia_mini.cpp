#include "../include/integral.h"

int test_integral_nuclear_ammonia_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/ammonia.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -48.24364223030854, -4.54235018816602, -0.00952939225634, -0.01186366108209, 0.00769117302004, -1.57041358951010, -1.57041416403629, -1.57041408155995, -4.54235018816602, -8.20160183703694, -0.12101335515777, -0.15065614649783, 0.09766987210038, -3.67308306137000, -3.67308374775032, -3.67308361112743, -0.00952939225634, -0.12101335515777, -8.14883504852827, 0.03625364222494, -0.02350326330237, -0.90277342756967, -2.19390737224789, 1.16156119524885, -0.01186366108209, -0.15065614649783, 0.03625364222494, -8.13282133126120, -0.02926021985611, -1.92815640004927, 0.85477173223541, -1.33575183363706, 0.00769117302004, 0.09766987210038, -0.02350326330237, -0.02926021985611, -8.15898525354944, -1.53419706697540, 1.15874216202089, 1.93728932875855, -1.57041358951010, -3.67308306137000, -0.90277342756967, -1.92815640004927, -1.53419706697540, -5.01669472298740, -1.80052481826866, -1.80052465718956, -1.57041416403629, -3.67308374775032, -2.19390737224789, 0.85477173223541, 1.15874216202089, -1.80052481826866, -5.01669473468874, -1.80052005241752, -1.57041408155995, -3.67308361112743, 1.16156119524885, -1.33575183363706, 1.93728932875855, -1.80052465718956, -1.80052005241752, -5.01669464767826;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}