#include "../include/integral.h"

int test_int_kinetic_ethane_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 15.89112181239595, -0.08588999412233, 0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000315060228, -0.00801237979553, 0.01080348303324, -0.00032310323026, -0.00029999653133, -0.00228505271246, -0.00228505428046, -0.00228506088898, -0.01038344653925, -0.01038344496713, -0.01038342780841, -0.08588999412233, 0.47224999256911, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00801237979553, 0.02125272243221, -0.10192489671585, 0.00304830055933, 0.00283030161451, -0.01500085308077, -0.01500085131464, -0.01500084387082, 0.10815381640004, 0.10815383044988, 0.10815398379498, 0.00000000000000, -0.00000000000000, 1.47772809059759, 0.00000000000000, 0.00000000000000, -0.01080348303324, 0.10192489671585, -0.22614061625249, 0.00805677889345, 0.00748059906366, -0.00115026662880, -0.00113956288159, -0.00110934432167, -0.09590691887673, -0.08350937444289, -0.10029677341098, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.47772809059759, -0.00000000000000, 0.00032310323026, -0.00304830055933, 0.00805677889345, 0.04300995871833, -0.00022372467419, 0.00001593604259, -0.00047565233103, 0.00056133399062, -0.20687257373552, 0.21976755594928, -0.00452948479827, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 1.47772809059759, 0.00029999653133, -0.00283030161451, 0.00748059906366, -0.00022372467419, 0.04304319036348, -0.00056729919285, 0.00034681863970, 0.00031483711113, 0.13226907043564, 0.11923894293526, -0.24374107751762, -0.00000315060228, -0.00801237979553, -0.01080348303324, 0.00032310323026, 0.00029999653133, 15.89112181239595, -0.08588999412233, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.01038345175442, -0.01038344761708, -0.01038341051474, -0.00228505639171, -0.00228505448281, -0.00228505700255, -0.00801237979553, 0.02125272243221, 0.10192489671585, -0.00304830055933, -0.00283030161451, -0.08588999412233, 0.47224999256911, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.10815376979264, 0.10815380676758, 0.10815413834537, -0.01500084893659, -0.01500085108673, -0.01500084824855, 0.01080348303324, -0.10192489671585, -0.22614061625249, 0.00805677889345, 0.00748059906366, -0.00000000000000, 0.00000000000000, 1.47772809059759, 0.00000000000000, 0.00000000000000, 0.10029901465933, 0.09590454505969, 0.08350957236215, 0.00113953934477, 0.00110943097045, 0.00115020059859, -0.00032310323026, 0.00304830055933, 0.00805677889345, 0.04300995871833, -0.00022372467419, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.47772809059759, 0.00000000000000, 0.00459685080844, 0.20683693747756, -0.21980012575470, 0.00047572667898, -0.00056130023139, -0.00001609897609, -0.00029999653133, 0.00283030161451, 0.00748059906366, -0.00022372467419, 0.04304319036348, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.47772809059759, 0.24373845869064, -0.13232647439269, -0.11918003681309, -0.00034667016731, -0.00031500573266, 0.00056727467546, -0.00228505271246, -0.01500085308077, -0.00115026662880, 0.00001593604259, -0.00056729919285, -0.01038345175442, 0.10815376979264, 0.10029901465933, 0.00459685080844, 0.24373845869064, 0.76003187992239, -0.00489386817738, -0.00489385576565, -0.01053751312109, -0.01053872446519, -0.00602063569977, -0.00228505428046, -0.01500085131464, -0.00113956288159, -0.00047565233103, 0.00034681863970, -0.01038344761708, 0.10815380676758, 0.09590454505969, 0.20683693747756, -0.13232647439269, -0.00489386817738, 0.76003187992239, -0.00489383254062, -0.00602063684628, -0.01053750863929, -0.01053873906450, -0.00228506088898, -0.01500084387082, -0.00110934432167, 0.00056133399062, 0.00031483711113, -0.01038341051474, 0.10815413834537, 0.08350957236215, -0.21980012575470, -0.11918003681309, -0.00489385576565, -0.00489383254062, 0.76003187992239, -0.01053875198854, -0.00602065115405, -0.01053753572638, -0.01038344653925, 0.10815381640004, -0.09590691887673, -0.20687257373552, 0.13226907043564, -0.00228505639171, -0.01500084893659, 0.00113953934477, 0.00047572667898, -0.00034667016731, -0.01053751312109, -0.00602063684628, -0.01053875198854, 0.76003187992239, -0.00489386135710, -0.00489387222048, -0.01038344496713, 0.10815383044988, -0.08350937444289, 0.21976755594928, 0.11923894293526, -0.00228505448281, -0.01500085108673, 0.00110943097045, -0.00056130023139, -0.00031500573266, -0.01053872446519, -0.01053750863929, -0.00602065115405, -0.00489386135710, 0.76003187992239, -0.00489384296771, -0.01038342780841, 0.10815398379498, -0.10029677341098, -0.00452948479827, -0.24374107751762, -0.00228505700255, -0.01500084824855, 0.00115020059859, -0.00001609897609, 0.00056727467546, -0.00602063569977, -0.01053873906450, -0.01053753572638, -0.00489387222048, -0.00489384296771, 0.76003187992239;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
