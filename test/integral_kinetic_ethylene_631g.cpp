#include "../include/integral.h"

int test_integral_kinetic_ethylene_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 16.20756797345014, -1.24756866124283, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.08995049652200, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00023071338674, -0.03244256917502, 0.02644267260909, 0.00036685617880, 0.00700418729300, 0.01132544314136, -0.06345309403249, -0.00088032552359, -0.01680758074247, -0.03490173568820, 0.02461470166613, -0.03490173071725, 0.02461461766223, -0.00045908869639, -0.00455284250212, -0.00045905437939, -0.00455285001897, -1.24756866124283, 0.93224556748206, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.30671845828079, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.03244256917502, 0.00420727678065, -0.13594864480064, -0.00188610285661, -0.03601034526616, 0.06988827388184, -0.24188028597938, -0.00335576054481, -0.06406972738837, 0.04608652299057, 0.11619451928785, 0.04608587345180, 0.11619423382416, -0.01727766228649, -0.01197215815282, -0.01727711199354, -0.01197248967178, -0.00000000000000, -0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, -0.02644267260909, 0.13594864480064, -0.35625805054525, -0.00556068835390, -0.10616722563139, 0.11155821609284, -0.15838054483492, -0.00377691747252, -0.07211064961323, -0.11582877746838, -0.06270509062618, -0.13718527772513, -0.07426700109662, -0.02409263592650, 0.01389425042111, -0.02287882080261, 0.01319384400332, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, -0.00036685617880, 0.00188610285661, -0.00556068835390, 0.04447435704604, -0.00147292610261, 0.00154771877542, -0.00377691747252, 0.11380396465592, -0.00100043735211, -0.20739866332641, -0.11227738273575, 0.20389724102644, 0.11038237392535, 0.01136019769001, -0.00655143887202, -0.01201210401735, 0.00692718509943, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, -0.00700418729300, 0.03601034526616, -0.10616722563139, -0.00147292610261, 0.01642971984985, 0.02954976038688, -0.07211064961323, -0.00100043735211, 0.09475555311123, -0.06308888828833, -0.03415381344852, -0.00402166879688, -0.00217718173481, -0.00454006974363, 0.00261826336229, -0.00789657229143, 0.00455382486156, 0.08995049652200, 0.30671845828079, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.25307171730000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.01132544314136, 0.06988827388184, -0.11155821609284, -0.00154771877542, -0.02954976038688, 0.09953708017092, -0.19641971885295, -0.00272505690193, -0.05202804267265, 0.09619926118246, 0.13492367142091, 0.09619898694262, 0.13492346269481, -0.01360421968194, 0.00952674917060, -0.01360439043723, 0.00952603544262, 0.00000000000000, -0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00000000000000, 0.06345309403249, 0.24188028597938, -0.15838054483492, -0.00377691747252, -0.07211064961323, 0.19641971885295, -0.08888467332693, -0.00401434842443, -0.07664379080702, -0.12793623211185, -0.09737323424612, -0.15152579646078, -0.11532753162937, 0.01859958159533, 0.07774090327890, 0.01766161273344, 0.07382458872898, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00088032552359, 0.00335576054481, -0.00377691747252, 0.11380396465592, -0.00100043735211, 0.00272505690193, -0.00401434842443, 0.20041035754875, -0.00106332853111, -0.22907781736933, -0.17435286003884, 0.22521142468795, 0.17141027013651, -0.00877010404835, -0.03665651332392, 0.00927290488870, 0.03876024234383, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.01680758074247, 0.06406972738837, -0.07211064961323, -0.00100043735211, 0.09475555311123, 0.05202804267265, -0.07664379080702, -0.00106332853111, 0.18016449265268, -0.06968350035411, -0.05303663935592, -0.00444206971517, -0.00338089584441, 0.00350494640365, 0.01464966821794, 0.00609586494585, 0.02548038672156, -0.00023071338674, -0.03244256917502, -0.02644267260909, -0.00036685617880, -0.00700418729300, 0.01132544314136, 0.06345309403249, 0.00088032552359, 0.01680758074247, 16.20756797345014, -1.24756866124283, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.08995049652200, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00045910114643, -0.00455283977442, -0.00045905240481, -0.00455285045140, -0.03490173012325, 0.02461460762535, -0.03490173548781, 0.02461469827928, -0.03244256917502, 0.00420727678065, 0.13594864480064, 0.00188610285661, 0.03601034526616, 0.06988827388184, 0.24188028597938, 0.00335576054481, 0.06406972738837, -1.24756866124283, 0.93224556748206, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.30671845828079, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.01727786192356, -0.01197203788110, -0.01727708032930, -0.01197250874734, 0.04608579584433, 0.11619419971664, 0.04608649680253, 0.11619450777857, 0.02644267260909, -0.13594864480064, -0.35625805054525, -0.00556068835390, -0.10616722563139, -0.11155821609284, -0.15838054483492, -0.00377691747252, -0.07211064961323, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.02287814815223, -0.01319399827626, 0.02409085201485, -0.01389277859231, 0.13720371107093, 0.07427702218286, 0.11585630798355, 0.06272000650721, 0.00036685617880, -0.00188610285661, -0.00556068835390, 0.04447435704604, -0.00147292610261, -0.00154771877542, -0.00377691747252, 0.11380396465592, -0.00100043735211, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.01201199290465, -0.00692740568965, -0.01136046400547, 0.00655138353084, -0.20388633293193, -0.11037653104104, 0.20740881893341, 0.11228290198063, 0.00700418729300, -0.03601034526616, -0.10616722563139, -0.00147292610261, 0.01642971984985, -0.02954976038688, -0.07211064961323, -0.00100043735211, 0.09475555311123, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00790166486683, -0.00455694892519, 0.00454484986078, -0.00262093647880, 0.00393512472797, 0.00213033120191, 0.06300468385130, 0.03410823501901, 0.01132544314136, 0.06988827388184, 0.11155821609284, 0.00154771877542, 0.02954976038688, 0.09953708017092, 0.19641971885295, 0.00272505690193, 0.05202804267265, 0.08995049652200, 0.30671845828079, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.25307171730000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.01360415773248, 0.00952700809892, -0.01360440026233, 0.00952599437418, 0.09619895417614, 0.13492343775597, 0.09619925012571, 0.13492366300553, -0.06345309403249, -0.24188028597938, -0.15838054483492, -0.00377691747252, -0.07211064961323, -0.19641971885295, -0.08888467332693, -0.00401434842443, -0.07664379080702, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00000000000000, -0.01766232024376, -0.07382192369383, -0.01859720194731, -0.07773555117028, 0.15154623948380, 0.11534310310854, 0.12796666394476, 0.09739639959954, -0.00088032552359, -0.00335576054481, -0.00377691747252, 0.11380396465592, -0.00100043735211, -0.00272505690193, -0.00401434842443, 0.20041035754875, -0.00106332853111, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, -0.00927346322072, -0.03875962414955, 0.00876983691547, 0.03665756323069, -0.22519949931968, -0.17140121165988, 0.22908907675012, 0.17436143582423, -0.01680758074247, -0.06406972738837, -0.07211064961323, -0.00100043735211, 0.09475555311123, -0.05202804267265, -0.07664379080702, -0.00106332853111, 0.18016449265268, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, -0.00610021993075, -0.02549664846001, -0.00350844754801, -0.01466516869955, 0.00434648122685, 0.00330814300649, 0.06959050694495, 0.05296586324756, -0.03490173568820, 0.04608652299057, -0.11582877746838, -0.20739866332641, -0.06308888828833, 0.09619926118246, -0.12793623211185, -0.22907781736933, -0.06968350035411, -0.00045910114643, -0.01727786192356, 0.02287814815223, 0.01201199290465, 0.00790166486683, -0.01360415773248, -0.01766232024376, -0.00927346322072, -0.00610021993075, 1.39567838062028, 0.25973499598884, -0.02394153676840, -0.00258626185921, -0.00232747729850, -0.01345635715637, -0.00008484051876, -0.00612021568080, 0.02461470166613, 0.11619451928785, -0.06270509062618, -0.11227738273575, -0.03415381344852, 0.13492367142091, -0.09737323424612, -0.17435286003884, -0.05303663935592, -0.00455283977442, -0.01197203788110, -0.01319399827626, -0.00692740568965, -0.00455694892519, 0.00952700809892, -0.07382192369383, -0.03875962414955, -0.02549664846001, 0.25973499598884, 0.24191663820000, -0.00258626185921, 0.03278171313617, -0.01345635715637, -0.00680024400675, -0.00612021568080, -0.01299270475463, -0.03490173071725, 0.04608587345180, -0.13718527772513, 0.20389724102644, -0.00402166879688, 0.09619898694262, -0.15152579646078, 0.22521142468795, -0.00444206971517, -0.00045905240481, -0.01727708032930, 0.02409085201485, -0.01136046400547, 0.00454484986078, -0.01360440026233, -0.01859720194731, 0.00876983691547, -0.00350844754801, -0.02394153676840, -0.00258626185921, 1.39567838062028, 0.25973499598884, -0.00008483840579, -0.00612017012884, -0.00232703513392, -0.01345604121700, 0.02461461766223, 0.11619423382416, -0.07426700109662, 0.11038237392535, -0.00217718173481, 0.13492346269481, -0.11532753162937, 0.17141027013651, -0.00338089584441, -0.00455285045140, -0.01197250874734, -0.01389277859231, 0.00655138353084, -0.00262093647880, 0.00952599437418, -0.07773555117028, 0.03665756323069, -0.01466516869955, -0.00258626185921, 0.03278171313617, 0.25973499598884, 0.24191663820000, -0.00612017012884, -0.01299268913489, -0.01345604121700, -0.00680149179080, -0.00045908869639, -0.01727766228649, -0.02409263592650, 0.01136019769001, -0.00454006974363, -0.01360421968194, 0.01859958159533, -0.00877010404835, 0.00350494640365, -0.03490173012325, 0.04608579584433, 0.13720371107093, -0.20388633293193, 0.00393512472797, 0.09619895417614, 0.15154623948380, -0.22519949931968, 0.00434648122685, -0.00232747729850, -0.01345635715637, -0.00008483840579, -0.00612017012884, 1.39567838062028, 0.25973499598884, -0.02394172474828, -0.00258609780418, -0.00455284250212, -0.01197215815282, 0.01389425042111, -0.00655143887202, 0.00261826336229, 0.00952674917060, 0.07774090327890, -0.03665651332392, 0.01464966821794, 0.02461460762535, 0.11619419971664, 0.07427702218286, -0.11037653104104, 0.00213033120191, 0.13492343775597, 0.11534310310854, -0.17140121165988, 0.00330814300649, -0.01345635715637, -0.00680024400675, -0.00612017012884, -0.01299268913489, 0.25973499598884, 0.24191663820000, -0.00258609780418, 0.03278200040993, -0.00045905437939, -0.01727711199354, -0.02287882080261, -0.01201210401735, -0.00789657229143, -0.01360439043723, 0.01766161273344, 0.00927290488870, 0.00609586494585, -0.03490173548781, 0.04608649680253, 0.11585630798355, 0.20740881893341, 0.06300468385130, 0.09619925012571, 0.12796666394476, 0.22908907675012, 0.06959050694495, -0.00008484051876, -0.00612021568080, -0.00232703513392, -0.01345604121700, -0.02394172474828, -0.00258609780418, 1.39567838062028, 0.25973499598884, -0.00455285001897, -0.01197248967178, 0.01319384400332, 0.00692718509943, 0.00455382486156, 0.00952603544262, 0.07382458872898, 0.03876024234383, 0.02548038672156, 0.02461469827928, 0.11619450777857, 0.06272000650721, 0.11228290198063, 0.03410823501901, 0.13492366300553, 0.09739639959954, 0.17436143582423, 0.05296586324756, -0.00612021568080, -0.01299270475463, -0.01345604121700, -0.00680149179080, -0.00258609780418, 0.03278200040993, 0.25973499598884, 0.24191663820000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}