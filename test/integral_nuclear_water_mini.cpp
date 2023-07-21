#include "../include/integral.h"

int test_integral_nuclear_water_mini(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "MINI", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -62.38414251822861, -5.93506627945969, -0.01520422781758, -0.00526742842217, -0.00608107563576, -1.76806043840923, -1.76805501262965, -5.93506627945969, -9.99815990336933, -0.19556835107280, -0.06775376624938, -0.07821939680467, -3.96960125199348, -3.96959241486472, -0.01520422781758, -0.19556835107280, -9.88811760828859, -0.00401214553216, -0.04049418314397, -1.91490595518672, -1.20995444039086, -0.00526742842217, -0.06775376624938, -0.00401214553216, -9.93529233823887, 0.03566183198654, 1.28978039605256, -2.37237634236364, -0.00608107563576, -0.07821939680467, -0.04049418314397, 0.03566183198654, -9.83534841371785, -1.32972443287039, 0.07990785502654, -1.76806043840923, -3.96960125199348, -1.91490595518672, 1.28978039605256, -1.32972443287039, -5.33074543603888, -2.15927818542771, -1.76805501262965, -3.96959241486472, -1.20995444039086, -2.37237634236364, 0.07990785502654, -2.15927818542771, -5.33074044832427;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}
