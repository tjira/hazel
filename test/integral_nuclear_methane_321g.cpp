#include "../include/integral.h"

int test_integral_nuclear_methane_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "3-21G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -36.39223078535572, -2.68326151108964, 0.00000011903156, 0.00000011226696, 0.00000020712638, -3.82593350325209, 0.00000001766028, 0.00000001665664, 0.00000003073058, -0.38152103106858, -1.71555380253914, -0.38152052443269, -1.71555326444447, -0.38151993926974, -1.71555264242405, -0.38151968517448, -1.71555237201956, -2.68326151108964, -8.19656974493183, 0.00000087345037, 0.00000082305303, 0.00000151640938, -5.78254163732090, 0.00000054344832, 0.00000050974650, 0.00000093273849, -1.37898541523534, -2.97378913167396, -1.37898468011967, -2.97378876900018, -1.37898360406288, -2.97378805082591, -1.37898300444166, -2.97378756445601, 0.00000011903156, 0.00000087345037, -8.71877986261150, 0.00000007079143, -0.00000056920872, 0.00000083856492, -3.49736232695250, 0.00000006420752, -0.00000047500967, -0.96090624457155, -0.74125394867006, -0.37277688412869, -0.28756410507390, -0.60786894816322, -0.46891725833772, 1.94156382203736, 1.49774718202755, 0.00000011226696, 0.00000082305303, 0.00000007079143, -8.71877952481203, 0.00000016335960, 0.00000078695325, 0.00000006420752, -3.49736203061622, 0.00000013608565, -0.29263214887335, -0.22573945419428, -1.49733514147738, -1.15506183267209, 1.69242515705719, 1.30555800636305, 0.09755338633772, 0.07525423012358, 0.00000020712638, 0.00000151640938, -0.00000056920872, 0.00000016335960, -8.71878076459993, 0.00000144105255, -0.00000047500967, 0.00000013608565, -3.49736307202793, -1.70063138416619, -1.31188634359884, 1.23294097992065, 0.95110610041571, 0.81688816958267, 0.63015827277342, -0.34917652385803, -0.26935850112337, -3.82593350325209, -5.78254163732090, 0.00000083856492, 0.00000078695325, 0.00000144105255, -6.05483796433210, 0.00000114598504, 0.00000105813738, 0.00000188998328, -1.94017308339983, -3.66836649021578, -1.94017287234526, -3.66836669950845, -1.94017223751434, -3.66836629520356, -1.94017173382037, -3.66836574293244, 0.00000001766028, 0.00000054344832, -3.49736232695250, 0.00000006420752, -0.00000047500967, 0.00000114598504, -4.45987213041489, 0.00000014471449, -0.00000088005934, -1.26645690647798, -1.20301307075711, -0.49131341954062, -0.46670027176992, -0.80116110087444, -0.76102632969781, 2.55894911659407, 2.43075935081981, 0.00000001665664, 0.00000050974650, 0.00000006420752, -3.49736203061622, 0.00000013608565, 0.00000105813738, 0.00000014471449, -4.45987150891079, 0.00000025092694, -0.38568385107473, -0.36636231340263, -1.97346148251475, -1.87460001862801, 2.23058796400963, 2.11884702387759, 0.12857373920781, 0.12213332280556, 0.00000003073058, 0.00000093273849, -0.00000047500967, 0.00000013608565, -3.49736307202793, 0.00000188998328, -0.00000088005934, 0.00000025092694, -4.45987347247378, -2.24140114921972, -2.12911709531048, 1.62499455542929, 1.54359132366171, 1.07664487621699, 1.02271127740235, -0.46020894330907, -0.43715376380802, -0.38152103106858, -1.37898541523534, -0.96090624457155, -0.29263214887335, -1.70063138416619, -1.94017308339983, -1.26645690647798, -0.38568385107473, -2.24140114921972, -5.54978216217835, -3.24944123427730, -0.05602198528300, -0.68548807131236, -0.05602307239011, -0.68549307268761, -0.05602263805970, -0.68549109836200, -1.71555380253914, -2.97378913167396, -0.74125394867006, -0.22573945419428, -1.31188634359884, -3.66836649021578, -1.20301307075711, -0.36636231340263, -2.12911709531048, -3.24944123427730, -4.27464548189909, -0.68548813845990, -1.86599597660813, -0.68549311183225, -1.86600330946260, -0.68549106379423, -1.86600032261332, -0.38152052443269, -1.37898468011967, -0.37277688412869, -1.49733514147738, 1.23294097992065, -1.94017287234526, -0.49131341954062, -1.97346148251475, 1.62499455542929, -0.05602198528300, -0.68548813845990, -5.54978309375851, -3.24944184937850, -0.05602408689933, -0.68549780359464, -0.05602399441631, -0.68549739342153, -1.71555326444447, -2.97378876900018, -0.28756410507390, -1.15506183267209, 0.95110610041571, -3.66836669950845, -0.46670027176992, -1.87460001862801, 1.54359132366171, -0.68548807131236, -1.86599597660813, -3.24944184937850, -4.27464658491591, -0.68549777559110, -1.86601034110647, -0.68549729170532, -1.86600963226201, -0.38151993926974, -1.37898360406288, -0.60786894816322, 1.69242515705719, 0.81688816958267, -1.94017223751434, -0.80116110087444, 2.23058796400963, 1.07664487621699, -0.05602307239011, -0.68549311183225, -0.05602408689933, -0.68549777559110, -5.54978308311887, -3.24944185776454, -0.05602395899395, -0.68549723056360, -1.71555264242405, -2.97378805082591, -0.46891725833772, 1.30555800636305, 0.63015827277342, -3.66836629520356, -0.76102632969781, 2.11884702387759, 1.02271127740235, -0.68549307268761, -1.86600330946260, -0.68549780359464, -1.86601034110647, -3.24944185776454, -4.27464681441492, -0.68549715685094, -1.86600948479989, -0.38151968517448, -1.37898300444166, 1.94156382203736, 0.09755338633772, -0.34917652385803, -1.94017173382037, 2.55894911659407, 0.12857373920781, -0.46020894330907, -0.05602263805970, -0.68549106379423, -0.05602399441631, -0.68549729170532, -0.05602395899395, -0.68549715685094, -5.54978244488834, -3.24944145203377, -1.71555237201956, -2.97378756445601, 1.49774718202755, 0.07525423012358, -0.26935850112337, -3.66836574293244, 2.43075935081981, 0.12213332280556, -0.43715376380802, -0.68549109836200, -1.86600032261332, -0.68549739342153, -1.86600963226201, -0.68549723056360, -1.86600948479989, -3.24944145203377, -4.27464630502537;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}