#include "../include/integral.h"

int test_int_kinetic_ethane_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix T = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Texp(system.shells.nbf(), system.shells.nbf()); Texp << 16.20756797345014, -1.24756866124283, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.08995049652200, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000784109811, -0.01869563954726, 0.02927140835536, -0.00087542939299, -0.00081282313739, 0.00191568013796, -0.03929603276211, 0.00117523904863, 0.00109119193204, -0.00025503359310, -0.00452747934787, -0.00025503482497, -0.00452748015884, -0.00025504001686, -0.00452748357664, -0.03488136936431, 0.02434991065124, -0.03488136962429, 0.02434991339912, -0.03488137246168, 0.02434994339045, -1.24756866124283, 0.93224556748206, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.30671845828079, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.01869563954726, -0.03967288292305, -0.01844397836651, 0.00055160997345, 0.00051216163500, 0.03088446351586, -0.18190023364365, 0.00544014859789, 0.00505109685220, -0.01338783235786, -0.01407906642823, -0.01338786088025, -0.01407905285767, -0.01338798109090, -0.01407899566279, 0.04404875641733, 0.11529386242571, 0.04404877746388, 0.11529387178105, 0.04404900717386, 0.11529397388861, -0.00000000000000, -0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, -0.02927140835536, 0.01844397836651, -0.15973198252030, 0.00470864383895, 0.00437190559137, 0.08129625321391, -0.19278962971526, 0.00780461438828, 0.00724646808077, -0.01993725802380, 0.00865010625537, -0.01975214877917, 0.00856981808732, -0.01923007418863, 0.00834340739835, -0.08788591830887, -0.04829673704427, -0.07652522612196, -0.04205358737452, -0.09190872388501, -0.05050733125937, -0.00000000000000, -0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00087542939299, -0.00055160997345, 0.00470864383895, -0.00243162547371, -0.00013075198199, -0.00243135310536, 0.00780461438828, 0.06793695566232, -0.00021672244383, 0.00027621508360, -0.00011984044243, -0.00824452582778, 0.00357703292691, 0.00973051745367, -0.00422180749363, -0.18957116262940, -0.10417674149471, 0.20138771276259, 0.11067037895868, -0.00415067358109, -0.00228095263048, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00081282313739, -0.00051216163500, 0.00437190559137, -0.00013075198199, -0.00241220381857, -0.00225747510313, 0.00724646808077, -0.00021672244383, 0.06796914720926, -0.00983284231801, 0.00426613984438, 0.00601143954533, -0.00260816906162, 0.00545758506720, -0.00236789807362, 0.12120698751712, 0.06660796310358, 0.10926661984403, 0.06004625634697, -0.22335645137064, -0.12274284530823, 0.08995049652200, 0.30671845828079, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.25307171730000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00191568013796, 0.03088446351586, -0.08129625321391, 0.00243135310536, 0.00225747510313, 0.06519260158249, -0.17643983083164, 0.00527684258059, 0.00489946965028, -0.01456170238424, 0.00444043906750, -0.01456169735129, 0.00444047679593, -0.01456167613902, 0.00444063580630, 0.09533438312127, 0.13426448962262, 0.09533439210114, 0.13426449647633, 0.09533449011078, 0.13426457128030, -0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00000000000000, 0.03929603276211, 0.18190023364365, -0.19278962971526, 0.00780461438828, 0.00724646808077, 0.17643983083164, -0.18015807314667, 0.00981290249254, 0.00911113363893, 0.00896469838157, 0.06720212148539, 0.00888151990398, 0.06657815141880, 0.00864699588486, 0.06481829894193, -0.09849012410492, -0.07517070171080, -0.08575865191971, -0.06545364739826, -0.10299813062980, -0.07861132419064, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, -0.00117523904863, -0.00544014859789, 0.00780461438828, 0.06793695566232, -0.00021672244383, -0.00527684258059, 0.00981290249254, 0.14765883518975, -0.00027248959442, -0.00012419886977, -0.00093103272185, 0.00370713693264, 0.02778964937309, -0.00437542484517, -0.03279839604282, -0.21244458376671, -0.16214426148343, 0.22568687000270, 0.17225117793835, -0.00465148031260, -0.00355015207154, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, -0.00109119193204, -0.00505109685220, 0.00724646808077, -0.00021672244383, 0.06796914720926, -0.00489946965028, 0.00911113363893, -0.00027248959442, 0.14769931028506, 0.00442129331472, 0.03314336721795, -0.00270303350639, -0.02026263252512, -0.00245405790713, -0.01839573664234, 0.13583167215696, 0.10367092338835, 0.12245057600632, 0.09345805520749, -0.25030591202707, -0.19104112935717, -0.00000784109811, -0.01869563954726, -0.02927140835536, 0.00087542939299, 0.00081282313739, 0.00191568013796, 0.03929603276211, -0.00117523904863, -0.00109119193204, 16.20756797345014, -1.24756866124283, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.08995049652200, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.03488136850189, 0.02434990153573, -0.03488136918607, 0.02434990876732, -0.03488137532126, 0.02434997361751, -0.00025503648363, -0.00452748125076, -0.00025503498394, -0.00452748026349, -0.00025503696352, -0.00452748156667, -0.01869563954726, -0.03967288292305, 0.01844397836651, -0.00055160997345, -0.00051216163500, 0.03088446351586, 0.18190023364365, -0.00544014859789, -0.00505109685220, -1.24756866124283, 0.93224556748206, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.30671845828079, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.04404868659978, 0.11529383139134, 0.04404874198798, 0.11529385601177, 0.04404923868952, 0.11529407679867, -0.01338789928435, -0.01407903458549, -0.01338786456098, -0.01407905110643, -0.01338791039554, -0.01407902929892, 0.02927140835536, -0.01844397836651, -0.15973198252030, 0.00470864383895, 0.00437190559137, -0.08129625321391, -0.19278962971526, 0.00780461438828, 0.00724646808077, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.09191066450521, 0.05050851815809, 0.08788373815024, 0.04829554414636, 0.07652554307686, 0.04205361726014, 0.01975230080300, -0.00856991705067, 0.01922992182889, -0.00834324379483, 0.01993726207951, -0.00865017570997, -0.00087542939299, 0.00055160997345, 0.00470864383895, -0.00243162547371, -0.00013075198199, 0.00243135310536, 0.00780461438828, 0.06793695566232, -0.00021672244383, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00421240042956, 0.00231487939654, 0.18953849623874, 0.10415880122082, -0.20141791552704, -0.11068659676682, 0.00824604828818, -0.00357770725195, -0.00972909523859, 0.00422114110504, -0.00027905524121, 0.00012107363888, -0.00081282313739, 0.00051216163500, 0.00437190559137, -0.00013075198199, -0.00241220381857, 0.00225747510313, 0.00724646808077, -0.00021672244383, 0.06796914720926, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.17833218887281, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.22335377650339, 0.12274166808532, -0.12125958387721, -0.06663687400620, -0.10921283372749, -0.06001649285728, -0.00600903641946, 0.00260713646389, -0.00546003832240, 0.00236893479122, 0.00983298382006, -0.00426623462427, 0.00191568013796, 0.03088446351586, 0.08129625321391, -0.00243135310536, -0.00225747510313, 0.06519260158249, 0.17643983083164, -0.00527684258059, -0.00489946965028, 0.08995049652200, 0.30671845828079, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.25307171730000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.09533435333241, 0.13426446688684, 0.09533437696474, 0.13426448492376, 0.09533458889072, 0.13426464667216, -0.01456169057460, 0.00444052759553, -0.01456169670180, 0.00444048166467, -0.01456168861393, 0.00444054229301, -0.03929603276211, -0.18190023364365, -0.19278962971526, 0.00780461438828, 0.00724646808077, -0.17643983083164, -0.18015807314667, 0.00981290249254, 0.00911113363893, 0.00000000000000, -0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00000000000000, 0.10300054290767, 0.07861320008667, 0.09848769111318, 0.07516884627297, 0.08575872260374, 0.06545365969315, -0.00888166239941, -0.06657862884099, -0.00864670836694, -0.06481788878270, -0.00896485226585, -0.06720206336659, 0.00117523904863, 0.00544014859789, 0.00780461438828, 0.06793695566232, -0.00021672244383, 0.00527684258059, 0.00981290249254, 0.14765883518975, -0.00027248959442, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.00000000000000, 0.00472066580657, 0.00360295815069, 0.21240799793594, 0.16211634127607, -0.22571996813872, -0.17227632982325, -0.00370785245503, -0.02779476648615, 0.00437467452810, 0.03279365453187, 0.00012547806221, 0.00094060498015, 0.00109119193204, 0.00505109685220, 0.00724646808077, -0.00021672244383, 0.06796914720926, 0.00489946965028, 0.00911113363893, -0.00027248959442, 0.14769931028506, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37667339168477, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42178619550000, 0.25030349159346, 0.19103936650769, -0.13589062883278, -0.10371592300738, -0.12238989409062, -0.09341168145315, 0.00270197549924, 0.02025452171132, 0.00245509885408, 0.01840403512191, -0.00442143193623, -0.03314380876987, -0.00025503359310, -0.01338783235786, -0.01993725802380, 0.00027621508360, -0.00983284231801, -0.01456170238424, 0.00896469838157, -0.00012419886977, 0.00442129331472, -0.03488136850189, 0.04404868659978, 0.09191066450521, 0.00421240042956, 0.22335377650339, 0.09533435333241, 0.10300054290767, 0.00472066580657, 0.25030349159346, 1.39567838062028, 0.25973499598884, -0.02883884222080, 0.00197061608907, -0.02883886792035, 0.00197064157418, -0.00173633476711, -0.01290963522738, -0.00173794343187, -0.01291152525780, -0.00008068674022, -0.00602893330891, -0.00452747934787, -0.01407906642823, 0.00865010625537, -0.00011984044243, 0.00426613984438, 0.00444043906750, 0.06720212148539, -0.00093103272185, 0.03314336721795, 0.02434990153573, 0.11529383139134, 0.05050851815809, 0.00231487939654, 0.12274166808532, 0.13426446688684, 0.07861320008667, 0.00360295815069, 0.19103936650769, 0.25973499598884, 0.24191663820000, 0.00197061608907, 0.04026236008256, 0.00197064157418, 0.04026239953624, -0.01290963522738, -0.00852533924537, -0.01291152525780, -0.00852048019978, -0.00602893330891, -0.01296038483115, -0.00025503482497, -0.01338786088025, -0.01975214877917, -0.00824452582778, 0.00601143954533, -0.01456169735129, 0.00888151990398, 0.00370713693264, -0.00270303350639, -0.03488136918607, 0.04404874198798, 0.08788373815024, 0.18953849623874, -0.12125958387721, 0.09533437696474, 0.09848769111318, 0.21240799793594, -0.13589062883278, -0.02883884222080, 0.00197061608907, 1.39567838062028, 0.25973499598884, -0.02883891600965, 0.00197068926229, -0.00008068681109, -0.00602893489668, -0.00173632881911, -0.01290962823416, -0.00173796283238, -0.01291154803555, -0.00452748015884, -0.01407905285767, 0.00856981808732, 0.00357703292691, -0.00260816906162, 0.00444047679593, 0.06657815141880, 0.02778964937309, -0.02026263252512, 0.02434990876732, 0.11529385601177, 0.04829554414636, 0.10415880122082, -0.06663687400620, 0.13426448492376, 0.07516884627297, 0.16211634127607, -0.10371592300738, 0.00197061608907, 0.04026236008256, 0.25973499598884, 0.24191663820000, 0.00197068926229, 0.04026247336246, -0.00602893489668, -0.01296038541110, -0.01290962823416, -0.00852535721334, -0.01291154803555, -0.00852042160544, -0.00025504001686, -0.01338798109090, -0.01923007418863, 0.00973051745367, 0.00545758506720, -0.01456167613902, 0.00864699588486, -0.00437542484517, -0.00245405790713, -0.03488137532126, 0.04404923868952, 0.07652554307686, -0.20141791552704, -0.10921283372749, 0.09533458889072, 0.08575872260374, -0.22571996813872, -0.12238989409062, -0.02883886792035, 0.00197064157418, -0.02883891600965, 0.00197068926229, 1.39567838062028, 0.25973499598884, -0.00173798000693, -0.01291156819952, -0.00008068769545, -0.00602895471101, -0.00173636476799, -0.01290967049967, -0.00452748357664, -0.01407899566279, 0.00834340739835, -0.00422180749363, -0.00236789807362, 0.00444063580630, 0.06481829894193, -0.03279839604282, -0.01839573664234, 0.02434997361751, 0.11529407679867, 0.04205361726014, -0.11068659676682, -0.06001649285728, 0.13426464667216, 0.06545365969315, -0.17227632982325, -0.09341168145315, 0.00197064157418, 0.04026239953624, 0.00197068926229, 0.04026247336246, 0.25973499598884, 0.24191663820000, -0.01291156819952, -0.00852036973420, -0.00602895471101, -0.01296039264846, -0.01290967049967, -0.00852524861762, -0.03488136936431, 0.04404875641733, -0.08788591830887, -0.18957116262940, 0.12120698751712, 0.09533438312127, -0.09849012410492, -0.21244458376671, 0.13583167215696, -0.00025503648363, -0.01338789928435, 0.01975230080300, 0.00824604828818, -0.00600903641946, -0.01456169057460, -0.00888166239941, -0.00370785245503, 0.00270197549924, -0.00173633476711, -0.01290963522738, -0.00008068681109, -0.00602893489668, -0.00173798000693, -0.01291156819952, 1.39567838062028, 0.25973499598884, -0.02883885634277, 0.00197063009321, -0.02883883384922, 0.00197060778735, 0.02434991065124, 0.11529386242571, -0.04829673704427, -0.10417674149471, 0.06660796310358, 0.13426448962262, -0.07517070171080, -0.16214426148343, 0.10367092338835, -0.00452748125076, -0.01407903458549, -0.00856991705067, -0.00357770725195, 0.00260713646389, 0.00444052759553, -0.06657862884099, -0.02779476648615, 0.02025452171132, -0.01290963522738, -0.00852533924537, -0.00602893489668, -0.01296038541110, -0.01291156819952, -0.00852036973420, 0.25973499598884, 0.24191663820000, 0.00197063009321, 0.04026238176246, 0.00197060778735, 0.04026234723060, -0.03488136962429, 0.04404877746388, -0.07652522612196, 0.20138771276259, 0.10926661984403, 0.09533439210114, -0.08575865191971, 0.22568687000270, 0.12245057600632, -0.00025503498394, -0.01338786456098, 0.01922992182889, -0.00972909523859, -0.00546003832240, -0.01456169670180, -0.00864670836694, 0.00437467452810, 0.00245509885408, -0.00173794343187, -0.01291152525780, -0.00173632881911, -0.01290962823416, -0.00008068769545, -0.00602895471101, -0.02883885634277, 0.00197063009321, 1.39567838062028, 0.25973499598884, -0.02883889441953, 0.00197066785228, 0.02434991339912, 0.11529387178105, -0.04205358737452, 0.11067037895868, 0.06004625634697, 0.13426449647633, -0.06545364739826, 0.17225117793835, 0.09345805520749, -0.00452748026349, -0.01407905110643, -0.00834324379483, 0.00422114110504, 0.00236893479122, 0.00444048166467, -0.06481788878270, 0.03279365453187, 0.01840403512191, -0.01291152525780, -0.00852048019978, -0.01290962823416, -0.00852535721334, -0.00602895471101, -0.01296039264846, 0.00197063009321, 0.04026238176246, 0.25973499598884, 0.24191663820000, 0.00197066785228, 0.04026244021752, -0.03488137246168, 0.04404900717386, -0.09190872388501, -0.00415067358109, -0.22335645137064, 0.09533449011078, -0.10299813062980, -0.00465148031260, -0.25030591202707, -0.00025503696352, -0.01338791039554, 0.01993726207951, -0.00027905524121, 0.00983298382006, -0.01456168861393, -0.00896485226585, 0.00012547806221, -0.00442143193623, -0.00008068674022, -0.00602893330891, -0.00173796283238, -0.01291154803555, -0.00173636476799, -0.01290967049967, -0.02883883384922, 0.00197060778735, -0.02883889441953, 0.00197066785228, 1.39567838062028, 0.25973499598884, 0.02434994339045, 0.11529397388861, -0.05050733125937, -0.00228095263048, -0.12274284530823, 0.13426457128030, -0.07861132419064, -0.00355015207154, -0.19104112935717, -0.00452748156667, -0.01407902929892, -0.00865017570997, 0.00012107363888, -0.00426623462427, 0.00444054229301, -0.06720206336659, 0.00094060498015, -0.03314380876987, -0.00602893330891, -0.01296038483115, -0.01291154803555, -0.00852042160544, -0.01290967049967, -0.00852524861762, 0.00197060778735, 0.04026234723060, 0.00197066785228, 0.04026244021752, 0.25973499598884, 0.24191663820000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC: " << T << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED KINETIC NORM: " << T.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC: " << Texp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED KINETIC NORM: " << Texp.norm() << std::endl;

    // return success or failure based on the error
    return (T - Texp).norm() > 1e-8;
}
