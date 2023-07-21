#include "../include/integral.h"

int test_int_overlap_ethylene_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix S = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Sexp(system.shells.nbf(), system.shells.nbf()); Sexp << 1.00000000000000, 0.24836239705696, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000424894366, 0.04686377296094, -0.07521947860024, -0.00104356813316, -0.01992428390226, 0.06345799591044, 0.06345780306567, 0.00634355959019, 0.00634339402542, 0.24836239705696, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.04686377296094, 0.41338471105430, -0.40702352996506, -0.00564689882491, -0.10781319568861, 0.49539895412283, 0.49539820228275, 0.10895921705033, 0.10895715907841, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.07521947860024, 0.40702352996506, -0.28831460561961, -0.00749586968758, -0.14311459980801, -0.22150376966432, -0.26234580234285, 0.12348262335127, 0.11726227526081, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00104356813316, 0.00564689882491, -0.00749586968758, 0.25187713247523, -0.00198552075246, -0.39661633968884, 0.38992219995904, -0.05822472131446, 0.06156640064172, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.01992428390226, 0.10781319568861, -0.14311459980801, -0.00198552075246, 0.21407264992213, -0.12064727682731, -0.00769082473551, 0.02326933938865, 0.04047280415560, 0.00000424894366, 0.04686377296094, 0.07521947860024, 0.00104356813316, 0.01992428390226, 1.00000000000000, 0.24836239705696, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00634361965439, 0.00634338449871, 0.06345778002439, 0.06345798813535, 0.04686377296094, 0.41338471105430, 0.40702352996506, 0.00564689882491, 0.10781319568861, 0.24836239705696, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.10895996364587, 0.10895704066098, 0.49539811245205, 0.49539892381028, -0.07521947860024, -0.40702352996506, -0.28831460561961, -0.00749586968758, -0.14311459980801, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.11725766851302, -0.12347442745428, 0.26238119504246, 0.22155645770454, -0.00104356813316, -0.00564689882491, -0.00749586968758, 0.25187713247523, -0.00198552075246, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.06156522253561, 0.05822653295226, -0.38990155054805, 0.39663583295000, -0.01992428390226, -0.10781319568861, -0.14311459980801, -0.00198552075246, 0.21407264992213, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.04049850510149, -0.02329401774911, 0.00752532654334, 0.12048627144989, 0.06345799591044, 0.49539895412283, -0.22150376966432, -0.39661633968884, -0.12064727682731, 0.00634361965439, 0.10895996364587, -0.11725766851302, -0.06156522253561, -0.04049850510149, 1.00000000000000, 0.15621473147726, 0.05352227853251, 0.01588447953226, 0.06345780306567, 0.49539820228275, -0.26234580234285, 0.38992219995904, -0.00769082473551, 0.00634338449871, 0.10895704066098, -0.12347442745428, 0.05822653295226, -0.02329401774911, 0.15621473147726, 1.00000000000000, 0.01588434408881, 0.05351820618095, 0.00634355959019, 0.10895921705033, 0.12348262335127, -0.05822472131446, 0.02326933938865, 0.06345778002439, 0.49539811245205, 0.26238119504246, -0.38990155054805, 0.00752532654334, 0.05352227853251, 0.01588434408881, 1.00000000000000, 0.15621543980837, 0.00634339402542, 0.10895715907841, 0.11726227526081, 0.06156640064172, 0.04047280415560, 0.06345798813535, 0.49539892381028, 0.22155645770454, 0.39663583295000, 0.12048627144989, 0.01588447953226, 0.05351820618095, 0.15621543980837, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP: " << S << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP NORM: " << S.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP: " << Sexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP NORM: " << Sexp.norm() << std::endl;

    // return success or failure based on the error
    return (S - Sexp).norm() > 1e-8;
}
