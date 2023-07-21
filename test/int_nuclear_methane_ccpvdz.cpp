#include "../include/integral.h"

int test_int_nuclear_methane_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "CC-PVDZ", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix V = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Vexp(system.shells.nbf(), system.shells.nbf()); Vexp << -35.94151538276792, 4.54872941599678, -3.53884195617617, 0.00000011200913, 0.00000010562535, 0.00000019482295, 0.00000002657214, 0.00000002504792, 0.00000004617338, 0.00000000403690, 0.00000000940789, -0.00000003560160, -0.00000003277441, -0.00000000966992, -1.48655521745671, -1.76477104422247, 0.94917786056303, 0.28906051233181, 1.67987422707136, -1.48655459978664, -1.76477068479980, 0.36822687403776, 1.47905864724803, -1.21789138026566, -1.48655387651335, -1.76477026263743, 0.60044908564780, -1.67176646501440, -0.80691663517608, -1.48655355669351, -1.76477007522261, -1.91786374412528, -0.09636250708840, 0.34491444528011, 4.54872941599678, -7.24950667509745, -5.22610294426887, 0.00000109794046, 0.00000102305663, 0.00000185328578, 0.00000093068843, 0.00000085995830, 0.00000153772369, 0.00000016567882, 0.00000033834160, -0.00000128295447, -0.00000118195097, -0.00000037700188, -3.07455705233169, -3.33809124332570, 1.09901831374913, 0.33469283197850, 1.94506483069416, -3.07455707690384, -3.33809144047867, 0.42635718703164, 1.71254987083896, -1.41015277018328, -3.07455656682327, -3.33809118421854, 0.69523968097120, -1.93568100483966, -0.93430074802571, -3.07455603134997, -3.33809079073446, -2.22062915761391, -0.11157469795119, 0.39936520067293, -3.53884195617617, -5.22610294426887, -5.58016475233140, 0.00000112108068, 0.00000104033292, 0.00000187270240, 0.00000107750609, 0.00000099131805, 0.00000176059272, 0.00000017403596, 0.00000034225756, -0.00000129861805, -0.00000119665681, -0.00000039055053, -3.34380034414249, -3.69768791553878, 0.98074452488808, 0.29867401161117, 1.73574147447397, -3.34380049765991, -3.69768822139134, 0.38047365764630, 1.52824920644788, -1.25839525927697, -3.34380007588717, -3.69768800673234, 0.62041959383820, -1.72736735031983, -0.83375333911335, -3.34379954340976, -3.69768758236282, -1.98164998408558, -0.09956724867394, 0.35638649095500, 0.00000011200913, 0.00000109794046, 0.00000112108068, -6.51857369141935, 0.00000011263855, -0.00000077966737, -3.99625533425486, 0.00000010921447, -0.00000069051630, 0.01277014032282, -0.01727650313387, -0.07549674996666, -0.05985840447933, 0.10908705590145, -1.25186280496026, -0.84365935093353, -0.60512699000400, 0.27968542737738, 1.62539003444656, -0.48565132515396, -0.32729139801079, -1.38530132455727, 0.55518125049041, -0.45714900385792, -0.79192855371500, -0.53369896173813, -1.15599336057434, -1.02326150254991, -0.49390068126728, 2.52946168828535, 1.70466450298728, 2.22595977518534, 0.18839127641921, -0.67431802472692, 0.00000010562535, 0.00000102305663, 0.00000104033292, 0.00000011263855, -6.51857318461979, 0.00000022302890, 0.00000010921447, -3.99625485879771, 0.00000019708557, -0.06549750825587, 0.12290477134375, -0.04255822082056, -0.01727650313387, -0.01180223774610, -0.38123906427149, -0.25692554778838, 0.27968517683735, -1.43834496990162, 0.49499241115231, -1.95072006675102, -1.31463569436822, 0.55518168678307, 0.70648056634032, -1.83623343673808, 2.20488422403486, 1.48592382599575, -1.02326205730913, 1.32544019212099, 1.37511546424230, 0.12709247578353, 0.08565095636872, 0.18839164542577, -1.51405256513556, -0.03388091739465, 0.00000019482295, 0.00000185328578, 0.00000187270240, -0.00000077966737, 0.00000022302890, -6.51857490356789, -0.00000069051630, 0.00000019708557, -3.99625639371884, -0.01727650313387, -0.04914062267203, -0.05459543246052, -0.08717465910301, -0.09138158791154, -2.21557228985869, -1.49312557006058, 1.62539010268545, 0.49499287534446, 1.35312860057945, 1.60626985850819, 1.08250368350866, -0.45714965983440, -1.83623462858646, -0.01152145440428, 1.06423890010763, 0.71721634318912, -0.49390116717600, 1.37511607159013, -0.85978744390206, -0.45490502916010, -0.30657060242379, -0.67431695108085, -0.03388079708635, -1.40224696668755, 0.00000002657214, 0.00000093068843, 0.00000107750609, -3.99625533425486, 0.00000010921447, -0.00000069051630, -3.99133656147198, 0.00000012850853, -0.00000075053368, 0.01180023240342, -0.01596441015645, -0.06976298050803, -0.05531246996232, 0.10080220602848, -1.23299806143064, -0.95462595379041, -0.85607103502745, 0.14324720595423, 0.83247972322615, -0.47833294894582, -0.37034012705020, -1.25565492709171, 0.28434856901933, -0.23413904813814, -0.77999497237464, -0.60389653243504, -1.13820973054490, -0.52408652983812, -0.25296235683474, 2.49134535460915, 1.92887924705301, 0.59393434104292, 0.09648877822073, -0.34536755462250, 0.00000002504792, 0.00000085995830, 0.00000099131805, 0.00000010921447, -3.99625485879771, 0.00000019708557, 0.00000012850853, -3.99133601710867, 0.00000021375835, -0.06052325427723, 0.11357043677170, -0.03932602633279, -0.01596441015645, -0.01090584066957, -0.37549403835452, -0.29071903254373, 0.14324699971922, -1.28282217558368, 0.25352136230221, -1.92132431299980, -1.48754991584802, 0.28434892801075, -0.18430174957779, -0.94046798759640, 2.17165878868982, 1.68136755178284, -0.52408697592458, 0.13271261486086, 0.70429638023690, 0.12517731951251, 0.09691657623218, 0.09648907155076, -1.32159786856664, -0.01735288474801, 0.00000004617338, 0.00000153772369, 0.00000176059272, -0.00000069051630, 0.00000019708557, -3.99625639371884, -0.00000075053368, 0.00000021375835, -3.99133769842446, -0.01596441015645, -0.04540860259520, -0.05044921127634, -0.08055407497610, -0.08444145336701, -2.18218509901274, -1.68951654837798, 0.83247978293044, 0.25352174548333, 0.14689320187512, 1.58206468630066, 1.22488537545061, -0.23413958504512, -0.94046895684766, -0.55204284888436, 1.04820185944197, 0.81155179403706, -0.25296274802882, 0.70429686992129, -0.98650128228172, -0.44805009415841, -0.34689400900655, -0.34536668622583, -0.01735278836230, -1.26433404872332, 0.00000000403690, 0.00000016567882, 0.00000017403596, 0.01277014032282, -0.06549750825587, -0.01727650313387, 0.01180023240342, -0.06052325427723, -0.01596441015645, -5.62627534790859, 0.05361824802059, -0.01998825339281, -0.01353376372162, -0.02705443644050, -0.19432597968458, -0.05313557785281, -0.10504924594387, -1.06033860640924, 0.42167075085611, -0.38574144071559, -0.10547544802282, -1.57703927096841, 0.29701456661288, -0.60856301989434, 0.71096548791852, 0.19440382125627, 1.43900153987190, 0.82013567791945, 0.74158000591776, -0.13089514376940, -0.03579153650573, -0.20953797016153, 2.26489760275236, 0.05770760363218, 0.00000000940789, 0.00000033834160, 0.00000034225756, -0.01727650313387, 0.12290477134375, -0.04914062267203, -0.01596441015645, 0.11357043677170, -0.04540860259520, 0.05361824802059, -5.72068147792175, 0.00898267173903, -0.02308019537691, -0.00834812847676, -0.34392228090261, -0.09404068959596, 0.42167088092430, -1.86681883132754, 0.40129146870543, 1.27582166039660, 0.34885625113873, -0.60856379172344, -0.98620337339588, 0.24708906807922, -0.95543516130904, -0.26125080444018, 0.74158056301760, -1.10071413352602, 0.99129283576103, 0.02354071463516, 0.00643707943350, 0.05770763248625, -0.40277045062340, 0.10223402047725, -0.00000003560160, -0.00000128295447, -0.00000129861805, -0.07549674996666, -0.04255822082056, -0.05459543246052, -0.06976298050803, -0.03932602633279, -0.05044921127634, -0.01998825339281, 0.00898267173903, -5.66774660802090, -0.07234697849470, 0.01982666545065, -0.95266494707235, -0.26049358950075, 1.81949327469661, 0.55330236580698, -0.24153053529064, -0.13153165156786, -0.03596612789773, 0.31274334410179, 1.26688435220340, 1.46519553328111, 0.37888516472938, 0.10360046732549, 0.11537838145978, -0.33159727520737, 1.50291802374419, 0.70529309314110, 0.19285255375742, 0.42847490535977, 0.01993942055425, -0.79128224420574, -0.00000003277441, -0.00000118195097, -0.00000119665681, -0.05985840447933, -0.01727650313387, -0.08717465910301, -0.05531246996232, -0.01596441015645, -0.08055407497610, -0.01353376372162, -0.02308019537691, -0.07234697849470, -5.67489489942852, -0.01185049164872, -1.12932648795717, -0.30879921484191, -0.61598933479156, 0.42167116628284, 1.32074831034097, 0.31762859176053, 0.08685087284019, 1.29654850910527, -0.60856344657106, 0.05844334930336, 0.34316371280389, 0.09383307585391, 0.69186093187956, 0.74158021161598, -0.36034007704342, 0.46851735823028, 0.12810914994670, 0.74942683785577, 0.05770744052062, 2.07123932779546, -0.00000000966992, -0.00000037700188, -0.00000039055053, 0.10908705590145, -0.01180223774610, -0.09138158791154, 0.10080220602848, -0.01090584066957, -0.08444145336701, -0.02705443644050, -0.00834812847676, 0.01982666545065, -0.01185049164872, -5.73925608137845, -0.28946140143778, -0.07914940930087, -0.77036774160860, 0.45200679387865, 0.62541480051516, 0.72668816060057, 0.19870274151250, -0.78025858341022, 0.37017866637680, 1.14137535788887, 0.86205417274916, 0.23571698474132, -1.38031726638192, -0.12663385524288, 0.89628729028065, -1.29928673325610, -0.35527256906343, -0.93579837245381, -0.27723546275507, 0.57609940300410, -1.48655521745671, -3.07455705233169, -3.34380034414249, -1.25186280496026, -0.38123906427149, -2.21557228985869, -1.23299806143064, -0.37549403835452, -2.18218509901274, -0.19432597968458, -0.34392228090261, -0.95266494707235, -1.12932648795717, -0.28946140143778, -4.67581880140837, -3.85755454568480, 0.51599847113242, 0.15714138621185, 0.91322445069043, -1.44062864424754, -1.88538776368243, -0.05250656517980, 0.44143602938917, -0.71856601667651, -1.44063480446723, -1.88539400584367, 0.02520434663515, -0.61296166473907, -0.58104240720731, -1.44063229869398, -1.88539149168475, -0.81752898274523, -0.08576364763328, -0.19558945233839, -1.76477104422247, -3.33809124332570, -3.69768791553878, -0.84365935093353, -0.25692554778838, -1.49312557006058, -0.95462595379041, -0.29071903254373, -1.68951654837798, -0.05313557785281, -0.09404068959596, -0.26049358950075, -0.30879921484191, -0.07914940930087, -3.85755454568480, -3.92118449867662, 0.44609269947271, 0.13585239721194, 0.78950381057776, -1.88538785659137, -2.31086064746004, -0.04049304770032, 0.58288067350594, -0.88575935058947, -1.88539409435413, -2.31086657593653, 0.06062623847013, -0.78912067644154, -0.70680811838501, -1.88539151314678, -2.31086415154129, -1.03595370798217, -0.10312198063019, -0.20525083528657, 0.94917786056303, 1.09901831374913, 0.98074452488808, -0.60512699000400, 0.27968517683735, 1.62539010268545, -0.85607103502745, 0.14324699971922, 0.83247978293044, -0.10504924594387, 0.42167088092430, 1.81949327469661, -0.61598933479156, -0.77036774160860, 0.51599847113242, 0.44609269947271, -4.65731698887787, -0.05334901825742, -0.31003749033274, 0.24945935172109, 0.33539784389671, -0.07685662630479, -0.11310790529347, 0.20027752645744, 0.20646690091126, 0.28626340203022, -0.10076229454159, 0.12392894551930, 0.12045592391788, 0.67270900052940, 0.81910948884480, 0.65598231021035, 0.09055509112937, 0.26840116896969, 0.28906051233181, 0.33469283197850, 0.29867401161117, 0.27968542737738, -1.43834496990162, 0.49499287534446, 0.14324720595423, -1.28282217558368, 0.25352174548333, -1.06033860640924, -1.86681883132754, 0.55330236580698, 0.42167116628284, 0.45200679387865, 0.15714138621185, 0.13585239721194, -0.05334901825742, -4.49838356740043, -0.09441826731757, -0.17710098246794, -0.18708188230481, -0.04233166378856, 0.01449467894624, -0.24600755987527, 0.40624513258507, 0.47959664524796, 0.02789918800392, 0.24606596312224, 0.38735771249767, 0.11457226522852, 0.14625819386976, 0.11580737430864, -0.08420269507990, 0.03806941339949, 1.67987422707136, 1.94506483069416, 1.73574147447397, 1.62539003444656, 0.49499241115231, 1.35312860057945, 0.83247972322615, 0.25352136230221, 0.14689320187512, 0.42167075085611, 0.40129146870543, -0.24153053529064, 1.32074831034097, 0.62541480051516, 0.91322445069043, 0.78950381057776, -0.31003749033274, -0.09441826731757, -5.03084735776503, 0.78763259673933, 0.98917533817252, 0.10347424104345, -0.40074906339257, 0.67664785577000, 0.71155037691450, 0.90222236909164, 0.02365063757261, 0.52783223173071, 0.46898773203526, 0.49829733858017, 0.65850640495305, 0.46200999505750, 0.05233865702606, 0.04250951721659, -1.48655459978664, -3.07455707690384, -3.34380049765991, -0.48565132515396, -1.95072006675102, 1.60626985850819, -0.47833294894582, -1.92132431299980, 1.58206468630066, -0.38574144071559, 1.27582166039660, -0.13153165156786, 0.31762859176053, 0.72668816060057, -1.44062864424754, -1.88538785659137, 0.24945935172109, -0.17710098246794, 0.78763259673933, -4.67581981297547, -3.85755549239509, 0.20017865747413, 0.80405695048625, -0.66207815863760, -1.44064068208166, -1.88540005593323, 0.13276286939172, -0.83328536726329, -0.04454676340114, -1.44064009659951, -1.88539947728453, -0.70997794052807, -0.30608395810927, 0.34090764652749, -1.76477068479980, -3.33809144047867, -3.69768822139134, -0.32729139801079, -1.31463569436822, 1.08250368350866, -0.37034012705020, -1.48754991584802, 1.22488537545061, -0.10547544802282, 0.34885625113873, -0.03596612789773, 0.08685087284019, 0.19870274151250, -1.88538776368243, -2.31086064746004, 0.33539784389671, -0.18708188230481, 0.98917533817252, -3.85755549239509, -3.92118556141730, 0.17305913070010, 0.69512600142285, -0.57238199795925, -1.88540005153436, -2.31087234374962, 0.18354940930419, -1.04091544542480, -0.09367387445177, -1.88539940583715, -2.31087174264046, -0.91303639713728, -0.35491405124928, 0.40788469882774, 0.36822687403776, 0.42635718703164, 0.38047365764630, -1.38530132455727, 0.55518168678307, -0.45714965983440, -1.25565492709171, 0.28434892801075, -0.23413958504512, -1.57703927096841, -0.60856379172344, 0.31274334410179, 1.29654850910527, -0.78025858341022, -0.05250656517980, -0.04049304770032, -0.07685662630479, -0.04233166378856, 0.10347424104345, 0.20017865747413, 0.17305913070010, -4.50850257024830, -0.10589899364387, 0.08719969300913, 0.01205631164212, 0.03329331150559, -0.09830869556159, -0.02546236317640, 0.00864703913166, 0.47830419505579, 0.56614340414426, 0.37075568226626, 0.26902939806160, -0.27782291608514, 1.47905864724803, 1.71254987083896, 1.52824920644788, 0.55518125049041, 0.70648056634032, -1.83623462858646, 0.28434856901933, -0.18430174957779, -0.94046895684766, 0.29701456661288, -0.98620337339588, 1.26688435220340, -0.60856344657106, 0.37017866637680, 0.44143602938917, 0.58288067350594, -0.11310790529347, 0.01449467894624, -0.40074906339257, 0.80405695048625, 0.69512600142285, -0.10589899364387, -4.90750255442663, 0.35025491651181, 0.80447503594809, 0.99777677513254, -0.10749088592827, 0.81910962722935, -0.03384744192363, 0.51279970510869, 0.66443664636230, 0.42183618980610, 0.12017718296835, -0.23095503342249, -1.21789138026566, -1.41015277018328, -1.25839525927697, -0.45714900385792, -1.83623343673808, -0.01152145440428, -0.23413904813814, -0.94046798759640, -0.55204284888436, -0.60856301989434, 0.24708906807922, 1.46519553328111, 0.05844334930336, 1.14137535788887, -0.71856601667651, -0.88575935058947, 0.20027752645744, -0.24600755987527, 0.67664785577000, -0.66207815863760, -0.57238199795925, 0.08719969300913, 0.35025491651181, -4.77054564016789, -0.25816572247141, -0.35958780389476, 0.03233075214908, -0.21002174267467, -0.10650436884892, -0.47141893124361, -0.60330367835010, -0.39831294694485, -0.20952771074760, 0.12336317598417, -1.48655387651335, -3.07455656682327, -3.34380007588717, -0.79192855371500, 2.20488422403486, 1.06423890010763, -0.77999497237464, 2.17165878868982, 1.04820185944197, 0.71096548791852, -0.95543516130904, 0.37888516472938, 0.34316371280389, 0.86205417274916, -1.44063480446723, -1.88539409435413, 0.20646690091126, 0.40624513258507, 0.71155037691450, -1.44064068208166, -1.88540005153436, 0.01205631164212, 0.80447503594809, -0.25816572247141, -4.67581993342120, -3.85755567994215, 0.32642118196379, -0.90881851677337, -0.43866173145230, -1.44063995113176, -1.88539935938493, -0.75297214041088, 0.27726837985267, 0.26481883216593, -1.76477026263743, -3.33809118421854, -3.69768800673234, -0.53369896173813, 1.48592382599575, 0.71721634318912, -0.60389653243504, 1.68136755178284, 0.81155179403706, 0.19440382125627, -0.26125080444018, 0.10360046732549, 0.09383307585391, 0.23571698474132, -1.88539400584367, -2.31086657593653, 0.28626340203022, 0.47959664524796, 0.90222236909164, -1.88540005593323, -2.31087234374962, 0.03329331150559, 0.99777677513254, -0.35958780389476, -3.85755567994215, -3.92118585130549, 0.28219874815454, -0.78569484138542, -0.37923330503813, -1.88539929233643, -2.31087168396702, -0.96217227369593, 0.31176858255974, 0.32092683674184, 0.60044908564780, 0.69523968097120, 0.62041959383820, -1.15599336057434, -1.02326205730913, -0.49390116717600, -1.13820973054490, -0.52408697592458, -0.25296274802882, 1.43900153987190, 0.74158056301760, 0.11537838145978, 0.69186093187956, -1.38031726638192, 0.02520434663515, 0.06062623847013, -0.10076229454159, 0.02789918800392, 0.02365063757261, 0.13276286939172, 0.18354940930419, -0.09830869556159, -0.10749088592827, 0.03233075214908, 0.32642118196379, 0.28219874815454, -4.55224221011203, 0.19518354839842, 0.09420976858481, 0.55601616013605, 0.66726334943640, 0.47777599748001, -0.29130378135940, -0.23499980734785, -1.67176646501440, -1.93568100483966, -1.72736735031983, -1.02326150254991, 1.32544019212099, 1.37511607159013, -0.52408652983812, 0.13271261486086, 0.70429686992129, 0.82013567791945, -1.10071413352602, -0.33159727520737, 0.74158021161598, -0.12663385524288, -0.61296166473907, -0.78912067644154, 0.12392894551930, 0.24606596312224, 0.52783223173071, -0.83328536726329, -1.04091544542480, -0.02546236317640, 0.81910962722935, -0.21002174267467, -0.90881851677337, -0.78569484138542, 0.19518354839842, -5.02556638549580, -0.26229806130999, -0.54160705045949, -0.70757195229688, -0.46936418719363, 0.11294491816547, 0.18060397084093, -0.80691663517608, -0.93430074802571, -0.83375333911335, -0.49390068126728, 1.37511546424230, -0.85978744390206, -0.25296235683474, 0.70429638023690, -0.98650128228172, 0.74158000591776, 0.99129283576103, 1.50291802374419, -0.36034007704342, 0.89628729028065, -0.58104240720731, -0.70680811838501, 0.12045592391788, 0.38735771249767, 0.46898773203526, -0.04454676340114, -0.09367387445177, 0.00864703913166, -0.03384744192363, -0.10650436884892, -0.43866173145230, -0.37923330503813, 0.09420976858481, -0.26229806130999, -4.60874220104630, -0.33388878996545, -0.42434762510215, -0.30812338909196, 0.14490692468769, 0.02357488221640, -1.48655355669351, -3.07455603134997, -3.34379954340976, 2.52946168828535, 0.12709247578353, -0.45490502916010, 2.49134535460915, 0.12517731951251, -0.44805009415841, -0.13089514376940, 0.02354071463516, 0.70529309314110, 0.46851735823028, -1.29928673325610, -1.44063229869398, -1.88539151314678, 0.67270900052940, 0.11457226522852, 0.49829733858017, -1.44064009659951, -1.88539940583715, 0.47830419505579, 0.51279970510869, -0.47141893124361, -1.44063995113176, -1.88539929233643, 0.55601616013605, -0.54160705045949, -0.33388878996545, -4.67581937467814, -3.85755523314358, -1.04260391253279, -0.05238507882641, 0.18750580465521, -1.76477007522261, -3.33809079073446, -3.69768758236282, 1.70466450298728, 0.08565095636872, -0.30657060242379, 1.92887924705301, 0.09691657623218, -0.34689400900655, -0.03579153650573, 0.00643707943350, 0.19285255375742, 0.12810914994670, -0.35527256906343, -1.88539149168475, -2.31086415154129, 0.81910948884480, 0.14625819386976, 0.65850640495305, -1.88539947728453, -2.31087174264046, 0.56614340414426, 0.66443664636230, -0.60330367835010, -1.88539935938493, -2.31087168396702, 0.66726334943640, -0.70757195229688, -0.42434762510215, -3.85755523314358, -3.92118543047684, -0.90135545082850, -0.04528811937773, 0.16210316719395, -1.91786374412528, -2.22062915761391, -1.98164998408558, 2.22595977518534, 0.18839164542577, -0.67431695108085, 0.59393434104292, 0.09648907155076, -0.34536668622583, -0.20953797016153, 0.05770763248625, 0.42847490535977, 0.74942683785577, -0.93579837245381, -0.81752898274523, -1.03595370798217, 0.65598231021035, 0.11580737430864, 0.46200999505750, -0.70997794052807, -0.91303639713728, 0.37075568226626, 0.42183618980610, -0.39831294694485, -0.75297214041088, -0.96217227369593, 0.47777599748001, -0.46936418719363, -0.30812338909196, -1.04260391253279, -0.90135545082850, -5.19733630379244, -0.03593486134566, 0.12862380208506, -0.09636250708840, -0.11157469795119, -0.09956724867394, 0.18839127641921, -1.51405256513556, -0.03388079708635, 0.09648877822073, -1.32159786856664, -0.01735278836230, 2.26489760275236, -0.40277045062340, 0.01993942055425, 0.05770744052062, -0.27723546275507, -0.08576364763328, -0.10312198063019, 0.09055509112937, -0.08420269507990, 0.05233865702606, -0.30608395810927, -0.35491405124928, 0.26902939806160, 0.12017718296835, -0.20952771074760, 0.27726837985267, 0.31176858255974, -0.29130378135940, 0.11294491816547, 0.14490692468769, -0.05238507882641, -0.04528811937773, -0.03593486134566, -4.48394299039682, 0.00646268822540, 0.34491444528011, 0.39936520067293, 0.35638649095500, -0.67431802472692, -0.03388091739465, -1.40224696668755, -0.34536755462250, -0.01735288474801, -1.26433404872332, 0.05770760363218, 0.10223402047725, -0.79128224420574, 2.07123932779546, 0.57609940300410, -0.19558945233839, -0.20525083528657, 0.26840116896969, 0.03806941339949, 0.04250951721659, 0.34090764652749, 0.40788469882774, -0.27782291608514, -0.23095503342249, 0.12336317598417, 0.26481883216593, 0.32092683674184, -0.23499980734785, 0.18060397084093, 0.02357488221640, 0.18750580465521, 0.16210316719395, 0.12862380208506, 0.00646268822540, -4.50526961485344;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR: " << V << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR NORM: " << V.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR: " << Vexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR NORM: " << Vexp.norm() << std::endl;

    // return success or failure based on the error
    return (V - Vexp).norm() > 1e-8;
}
