#include "../include/integral.h"

int test_int_kinetic_methane_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "3-21G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 16.57897377972110, -1.43502704208436, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.10274966905493, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.04456009635039, 0.02191642501858, -0.04456007523532, 0.02191640875950, -0.04456005123644, 0.02191639027985, -0.04456004104223, 0.02191638243011, -1.43502704208436, 1.35736735247949, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.34641224940510, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.01250517623602, 0.10494538494962, -0.01250527928300, 0.10494532765394, -0.01250539640343, 0.10494526253314, -0.01250544615356, 0.10494523487126, -0.00000000000000, 0.00000000000000, 2.58667015032044, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.00000000000000, 0.00000000000000, 0.08206309474049, 0.06404677686978, 0.03183577698561, 0.02484651387296, 0.05191296492923, 0.04051599874858, -0.16581250153628, -0.12941013393816, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.58667015032044, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.00000000000000, 0.02499131138049, 0.01950466222035, 0.12787492133676, 0.09980111396161, -0.14453574476935, -0.11280438445740, -0.00833120171169, -0.00650217516404, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.58667015032044, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.14523692947742, 0.11335128470243, -0.10529519470746, -0.08217856649885, -0.06976351140789, -0.05444763836456, 0.02982022243953, 0.02327351040610, 0.10274966905493, 0.34641224940510, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.29378550000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.09055667970820, 0.14041509785259, 0.09055662252349, 0.14041505244862, 0.09055655752881, 0.14041500084364, 0.09055652992051, 0.14041497892297, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.48964250000000, 0.00000000000000, 0.00000000000000, 0.14357957613541, 0.11400041565847, 0.05570072354596, 0.04422569725812, 0.09082845071385, 0.07211670140465, -0.29011062661160, -0.23034438246660, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.48964250000000, 0.00000000000000, 0.04372540307463, 0.03471743168151, 0.22373336906647, 0.17764157477609, -0.25288399127402, -0.20078685858223, -0.01457652545262, -0.01157358761073, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.41750973612682, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.48964250000000, 0.25411004592893, 0.20176024779158, -0.18422727781289, -0.14627421865565, -0.12206029199400, -0.09691440910765, 0.05217437369002, 0.04142583134088, -0.04456009635039, -0.01250517623602, 0.08206309474049, 0.02499131138049, 0.14523692947742, 0.09055667970820, 0.14357957613541, 0.04372540307463, 0.25411004592893, 1.54940494556963, 0.29315085381580, -0.02145068937070, -0.00677520644224, -0.02145101959973, -0.00677496973990, -0.02145089078418, -0.00677506207276, 0.02191642501858, 0.10494538494962, 0.06404677686978, 0.01950466222035, 0.11335128470243, 0.14041509785259, 0.11400041565847, 0.03471743168151, 0.20176024779158, 0.29315085381580, 0.27478737000000, -0.00677520644224, 0.03140748415976, -0.00677496973990, 0.03140796208836, -0.00677506207276, 0.03140777565850, -0.04456007523532, -0.01250527928300, 0.03183577698561, 0.12787492133676, -0.10529519470746, 0.09055662252349, 0.05570072354596, 0.22373336906647, -0.18422727781289, -0.02145068937070, -0.00677520644224, 1.54940494556963, 0.29315085381580, -0.02145132585193, -0.00677475022254, -0.02145130015054, -0.00677476864499, 0.02191640875950, 0.10494532765394, 0.02484651387296, 0.09980111396161, -0.08217856649885, 0.14041505244862, 0.04422569725812, 0.17764157477609, -0.14627421865565, -0.00677520644224, 0.03140748415976, 0.29315085381580, 0.27478737000000, -0.00677475022254, 0.03140840531432, -0.00677476864499, 0.03140836811784, -0.04456005123644, -0.01250539640343, 0.05191296492923, -0.14453574476935, -0.06976351140789, 0.09055655752881, 0.09082845071385, -0.25288399127402, -0.12206029199400, -0.02145101959973, -0.00677496973990, -0.02145132585193, -0.00677475022254, 1.54940494556963, 0.29315085381580, -0.02145129178221, -0.00677477464330, 0.02191639027985, 0.10494526253314, 0.04051599874858, -0.11280438445740, -0.05444763836456, 0.14041500084364, 0.07211670140465, -0.20078685858223, -0.09691440910765, -0.00677496973990, 0.03140796208836, -0.00677475022254, 0.03140840531432, 0.29315085381580, 0.27478737000000, -0.00677477464330, 0.03140835600674, -0.04456004104223, -0.01250544615356, -0.16581250153628, -0.00833120171169, 0.02982022243953, 0.09055652992051, -0.29011062661160, -0.01457652545262, 0.05217437369002, -0.02145089078418, -0.00677506207276, -0.02145130015054, -0.00677476864499, -0.02145129178221, -0.00677477464330, 1.54940494556963, 0.29315085381580, 0.02191638243011, 0.10494523487126, -0.12941013393816, -0.00650217516404, 0.02327351040610, 0.14041497892297, -0.23034438246660, -0.01157358761073, 0.04142583134088, -0.00677506207276, 0.03140777565850, -0.00677476864499, 0.03140836811784, -0.00677477464330, 0.03140835600674, 0.29315085381580, 0.27478737000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
