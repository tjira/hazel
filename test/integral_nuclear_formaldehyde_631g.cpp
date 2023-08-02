#include "../include/integral.h"

int test_integral_nuclear_formaldehyde_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -64.61071832929831, -5.65517589772133, -0.09064364796294, -0.00229319450161, -0.01043675108548, -5.84190580994243, -0.01249053263203, -0.00031599810243, -0.00143816564023, -0.00003744068900, -0.68958432337121, 1.46254666979001, 0.03700097041781, 0.16839835221428, -1.72040572126900, 3.18689335293939, 0.08062521977389, 0.36694048838183, -0.00105212255578, -0.37995830806237, -0.00105212838453, -0.37995886133809, -5.65517589772133, -13.66067904943506, -0.61841001758019, -0.01564516084639, -0.07120401227395, -9.41937125330113, -0.40892452794348, -0.01034538548578, -0.04708375712839, -0.14473487656499, -2.11016479641508, 3.43368615389889, 0.08686883012348, 0.39535634010704, -3.49551162789840, 5.58615475526203, 0.14132413636978, 0.64319264811493, -0.02881031686374, -0.92574571671630, -0.02881041527748, -0.92574687848890, -0.09064364796294, -0.61841001758019, -14.94116053575950, -0.00682840255451, -0.03033759078996, -0.59995120655179, -5.52280141367957, -0.00691216680504, -0.03069748402368, -0.36431436345620, -2.44481346660678, 3.41421877843824, 0.11067563214871, 0.50401854435692, -1.83499963185664, 1.25319744217390, 0.07179097139476, 0.32711434207792, -0.04904270842088, -0.63044690979507, -0.05396793095627, -0.69181919713562, -0.00229319450161, -0.01564516084639, -0.00682840255451, -14.67447240634586, -0.00009800258554, -0.01517817122562, -0.00691216680504, -5.25289230804073, -0.00008779962897, -0.00921678966785, -0.06185134103581, 0.11067563573270, -0.95897406852773, 0.01303443302286, -0.04642366079055, 0.07179097571733, -1.58425208855617, 0.00861973803784, 0.00289804501799, 0.03562404414908, -0.00550412958242, -0.06907616682613, -0.01043675108548, -0.07120401227395, -0.03033759078996, -0.00009800258554, -14.68132155119667, -0.06907865626254, -0.03069748402368, -0.00008779962897, -5.25988253915101, -0.04194733713943, -0.28149703314477, 0.50401852582282, 0.01303443214130, -0.90523413580860, -0.21128275980468, 0.32711431968319, 0.00861973697358, -1.55021773720920, -0.02824060792979, -0.35413371742076, 0.01638000914059, 0.20188791635146, -5.84190580994243, -9.41937125330113, -0.59995120655179, -0.01517817122562, -0.06907865626254, -9.71900808660974, -1.06214468225704, -0.02687120934018, -0.12229581913063, -1.44419756019607, -3.91628433011064, 3.80466881058730, 0.09625431455930, 0.43807150048216, -5.24493255723057, 6.04066857795453, 0.15282287481853, 0.69552563445037, -0.38166026540360, -1.90851689363436, -0.38166090512503, -1.90851858807824, -0.01249053263203, -0.40892452794348, -5.52280141367957, -0.00691216680504, -0.03069748402368, -1.06214468225704, -8.05081061241346, -0.02505056337506, -0.11007992791002, -3.29373668662097, -5.83845897771221, 3.10904686498738, 0.14015475665636, 0.64289486324343, -4.80471399231863, 0.96711691758666, 0.11727803204725, 0.54057378270999, -0.89848141407129, -2.36097961782779, -0.99138671747658, -2.57120428273457, -0.00031599810243, -0.01034538548578, -0.00691216680504, -5.25289230804073, -0.00008779962897, -0.02687120934018, -0.02505056337506, -7.07745204139658, 0.00077181827575, -0.08332824968323, -0.14770718935056, 0.14015477441861, -2.44803529760481, 0.02081177249309, -0.12155447769468, 0.11727805281644, -3.69368449280600, 0.01984790620791, 0.05534124400136, 0.11693092078585, -0.10315322470075, -0.24171056490055, -0.00143816564023, -0.04708375712839, -0.03069748402368, -0.00008779962897, -5.25988253915101, -0.12229581913063, -0.11007992791002, 0.00077181827575, -7.10823958830387, -0.37924248181117, -0.67224304049401, 0.64289477239996, 0.02081176814968, -2.40152486526459, -0.55321710550786, 0.54057367559103, 0.01984790110653, -3.66694028198006, -0.52965048466261, -1.23624826042611, 0.31205120144851, 0.66835578576803, -0.00003744068900, -0.14473487656499, -0.36431436345620, -0.00921678966785, -0.04194733713943, -1.44419756019607, -3.29373668662097, -0.08332824968323, -0.37924248181117, -38.59305990970996, -3.38555861569718, 0.11900328172155, 0.00301066150394, 0.01370211025570, -4.07106242610478, 0.01999183171669, 0.00050577292720, 0.00230187166633, -0.63612646359122, -1.98521538170805, -0.63612644738161, -1.98521536916219, -0.68958432337121, -2.11016479641508, -2.44481346660678, -0.06185134103581, -0.28149703314477, -3.91628433011064, -5.83845897771221, -0.14770718935056, -0.67224304049401, -3.38555861569718, -10.15593952482487, 0.85552014342860, 0.02164378624859, 0.09850510855273, -7.71265971305316, 0.59954177285671, 0.01516779535720, 0.06903160091189, -2.22192643027363, -4.20376785257915, -2.22192656789909, -4.20376805101338, 1.46254666979001, 3.43368615389889, 3.41421877843824, 0.11067563573270, 0.50401852582282, 3.80466881058730, 3.10904686498738, 0.14015477441861, 0.64289477239996, 0.11900328172155, 0.85552014342860, -11.07988477924484, -0.02020119541587, -0.06354555127049, 0.83447283589656, -5.23850610521795, -0.01860285709485, -0.06029614785319, -1.23009994715514, -0.59118032676472, -1.75835592922333, -0.97443762343393, 0.03700097041781, 0.08686883012348, 0.11067563214871, -0.95897406852773, 0.01303443214130, 0.09625431455930, 0.14015475665636, -2.44803529760481, 0.02081176814968, 0.00301066150394, 0.02164378624859, -0.02020119541587, -10.39885554092802, 0.02409065580256, 0.02111131147240, -0.01860285709485, -4.60403521611573, 0.02053003635102, 0.41280158243891, 0.30711554026559, -0.48840692659506, -0.34672440537368, 0.16839835221428, 0.39535634010704, 0.50401854435692, 0.01303443302286, -0.90523413580860, 0.43807150048216, 0.64289486324343, 0.02081177249309, -2.40152486526459, 0.01370211025570, 0.09850510855273, -0.06354555127049, 0.02409065580256, -10.54111046890088, 0.09608170582815, -0.06029614785319, 0.02053003635102, -4.72675601050454, -2.56502798501253, -1.82627631355546, 2.22093762018872, 1.64601172145414, -1.72040572126900, -3.49551162789840, -1.83499963185664, -0.04642366079055, -0.21128275980468, -5.24493255723057, -4.80471399231863, -0.12155447769468, -0.55321710550786, -4.07106242610478, -7.71265971305316, 0.83447283589656, 0.02111131147240, 0.09608170582815, -8.08419705680052, 1.07845012326693, 0.02728369031712, 0.12417338654602, -2.65841424093863, -4.94772563851828, -2.65841446301488, -4.94772602615583, 3.18689335293939, 5.58615475526203, 1.25319744217390, 0.07179097571733, 0.32711431968319, 6.04066857795453, 0.96711691758666, 0.11727805281644, 0.54057367559103, 0.01999183171669, 0.59954177285671, -5.23850610521795, -0.01860285709485, -0.06029614785319, 1.07845012326693, -7.10691631300662, -0.03430330625044, -0.11677879041180, -1.42641339192879, -1.01017100772038, -2.03175157527765, -1.60334672659583, 0.08062521977389, 0.14132413636978, 0.07179097139476, -1.58425208855617, 0.00861973697358, 0.15282287481853, 0.11727803204725, -3.69368449280600, 0.01984790110653, 0.00050577292720, 0.01516779535720, -0.01860285709485, -4.60403521611573, 0.02053003635102, 0.02728369031712, -0.03430330625044, -5.91392260456540, 0.03265265439844, 0.47261130636860, 0.47292128321794, -0.56009989427178, -0.53904113627911, 0.36694048838183, 0.64319264811493, 0.32711434207792, 0.00861973803784, -1.55021773720920, 0.69552563445037, 0.54057378270999, 0.01984790620791, -3.66694028198006, 0.00230187166633, 0.06903160091189, -0.06029614785319, 0.02053003635102, -4.72675601050454, 0.12417338654602, -0.11677879041180, 0.03265265439844, -6.11417660701308, -2.94124872282969, -2.83752789012994, 2.54307618490103, 2.53660843313942, -0.00105212255578, -0.02881031686374, -0.04904270842088, 0.00289804501799, -0.02824060792979, -0.38166026540360, -0.89848141407129, 0.05534124400136, -0.52965048466261, -0.63612646359122, -2.22192643027363, -1.23009994715514, 0.41280158243891, -2.56502798501253, -2.65841424093863, -1.42641339192879, 0.47261130636860, -2.94124872282969, -6.89987672744151, -4.18304908550170, -0.11365574859280, -0.95070802886904, -0.37995830806237, -0.92574571671630, -0.63044690979507, 0.03562404414908, -0.35413371742076, -1.90851689363436, -2.36097961782779, 0.11693092078585, -1.23624826042611, -1.98521538170805, -4.20376785257915, -0.59118032676472, 0.30711554026559, -1.82627631355546, -4.94772563851828, -1.01017100772038, 0.47292128321794, -2.83752789012994, -4.18304908550170, -5.61089940365182, -0.95070808991821, -2.47876861470137, -0.00105212838453, -0.02881041527748, -0.05396793095627, -0.00550412958242, 0.01638000914059, -0.38166090512503, -0.99138671747658, -0.10315322470075, 0.31205120144851, -0.63612644738161, -2.22192656789909, -1.75835592922333, -0.48840692659506, 2.22093762018872, -2.65841446301488, -2.03175157527765, -0.56009989427178, 2.54307618490103, -0.11365574859280, -0.95070808991821, -6.89987734736981, -4.18304949445143, -0.37995886133809, -0.92574687848890, -0.69181919713562, -0.06907616682613, 0.20188791635146, -1.90851858807824, -2.57120428273457, -0.24171056490055, 0.66835578576803, -1.98521536916219, -4.20376805101338, -0.97443762343393, -0.34672440537368, 1.64601172145414, -4.94772602615583, -1.60334672659583, -0.53904113627911, 2.53660843313942, -0.95070802886904, -2.47876861470137, -4.18304949445143, -5.61090001835214;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}