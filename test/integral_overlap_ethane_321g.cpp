#include "../include/integral.h"

int test_integral_overlap_ethane_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethane.xyz", "3-21G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.19144744450067, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.18031400188717, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000003514, 0.00165334597581, -0.00528984891536, 0.00015820520724, 0.00014689117584, 0.03627126005389, -0.09038250777831, 0.00270309863333, 0.00250978677370, 0.00000201298284, 0.00852929429132, 0.00000201299747, 0.00852930957856, 0.00000201305911, 0.00852937400828, 0.01936364917292, 0.08126592813781, 0.01936365384454, 0.08126593293284, 0.01936370483231, 0.08126598526753, 0.19144744450067, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.76135740608460, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00165334597581, 0.05420652572519, -0.10040316485566, 0.00300278963678, 0.00278804540156, 0.21052066164641, -0.42200122290709, 0.01262092585105, 0.01171834145538, 0.00132075867257, 0.06414574506885, 0.00132076396622, 0.06414583914144, 0.00132078627698, 0.06414623562007, 0.19687100392480, 0.40307559767819, 0.19687102915351, 0.40307561699107, 0.19687130450898, 0.40307582777895, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.52895881204315, -0.00000000000000, 0.00000000000000, 0.00528984891536, 0.10040316485566, -0.16650093302452, 0.00586404765501, 0.00544468080575, 0.17951793448373, -0.23126535356657, 0.01109499526762, 0.01030153766262, 0.00313884064795, 0.06382323888776, 0.00310970417580, 0.06323063664233, 0.00302753686096, 0.06155924802670, -0.10777995984544, -0.08553323626446, -0.09384762989784, -0.07447665397417, -0.11271326806701, -0.08944817024160, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.52895881204315, 0.00000000000000, -0.00015820520724, -0.00300278963678, 0.00586404765501, 0.02939767842756, -0.00016283581422, -0.00536890041324, 0.01109499526762, 0.13938207204499, -0.00030809138916, -0.00004348617703, -0.00088422095174, 0.00129798720538, 0.02639240027670, -0.00153194938190, -0.03114929934657, -0.23248289019685, -0.18449639434980, 0.24697423962126, 0.19599658515435, -0.00509022391160, -0.00403955295434, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.52895881204315, -0.00014689117584, -0.00278804540156, 0.00544468080575, -0.00016283581422, 0.02942186575520, -0.00498494397517, 0.01030153766262, -0.00030809138916, 0.13942783536437, 0.00154804262029, 0.03147693847665, -0.00094641848164, -0.01924383791544, -0.00085922913249, -0.01747080273757, 0.14864365644642, 0.11796230954508, 0.13400043121686, 0.10634156407564, -0.27391562534731, -0.21737681737109, 0.18031400188717, 0.76135740608460, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.03627126005389, 0.21052066164641, -0.17951793448373, 0.00536890041324, 0.00498494397517, 0.43742292602199, -0.56204331679435, 0.01680920964511, 0.01560710050444, 0.04517024141366, 0.20221517434217, 0.04517031385342, 0.20221536711673, 0.04517061915832, 0.20221617958385, 0.34031696524051, 0.67062738472784, 0.34031698325398, 0.67062740577061, 0.34031717985950, 0.67062763543913, 0.00000000000000, -0.00000000000000, 0.52895881204315, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.09038250777831, 0.42200122290709, -0.23126535356657, 0.01109499526762, 0.01030153766262, 0.56204331679435, -0.28474473177771, 0.02159809964134, 0.02005351345626, 0.11975617197283, 0.31869173496317, 0.11864424188648, 0.31573249216470, 0.11550814510939, 0.30738596067959, -0.18336915941501, -0.21420623256681, -0.15966567009192, -0.18651653871154, -0.19176214056360, -0.22401059818092, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.52895881204315, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00270309863333, -0.01262092585105, 0.01109499526762, 0.13938207204499, -0.00030809138916, -0.01680920964511, 0.02159809964134, 0.43677698473249, -0.00059974685533, -0.00165912790081, -0.00441522420535, 0.04952197998742, 0.13178640538930, -0.05844772157416, -0.15553889319746, -0.39552985745125, -0.46204585821634, 0.42018437234376, 0.49084649633901, -0.00866013602460, -0.01011650290040, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.52895881204315, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00250978677370, -0.01171834145538, 0.01030153766262, -0.00030809138916, 0.13942783536437, -0.01560710050444, 0.02005351345626, -0.00059974685533, 0.43686607000996, 0.05906246256191, 0.15717535351181, -0.03610861256037, -0.09609115496058, -0.03278175225457, -0.08723757445839, 0.25289174698190, 0.29542038879956, 0.22797878503845, 0.26631782436734, -0.46602008398150, -0.54439024027476, 0.00000000003514, 0.00165334597581, 0.00528984891536, -0.00015820520724, -0.00014689117584, 0.03627126005389, 0.09038250777831, -0.00270309863333, -0.00250978677370, 1.00000000000000, 0.19144744450067, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.18031400188717, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.01936363367582, 0.08126591223131, 0.01936364597010, 0.08126592485038, 0.01936375622093, 0.08126603801354, 0.00000201301716, 0.00852933016215, 0.00000201299935, 0.00852931155134, 0.00000201302286, 0.00852933611745, 0.00165334597581, 0.05420652572519, 0.10040316485566, -0.00300278963678, -0.00278804540156, 0.21052066164641, 0.42200122290709, -0.01262092585105, -0.01171834145538, 0.19144744450067, 1.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.76135740608460, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.19687092023383, 0.40307553361171, 0.19687098662820, 0.40307558443743, 0.19687158202877, 0.40307604022348, 0.00132077109389, 0.06414596580588, 0.00132076464935, 0.06414585128124, 0.00132077315609, 0.06414600245280, -0.00528984891536, -0.10040316485566, -0.16650093302452, 0.00586404765501, 0.00544468080575, -0.17951793448373, -0.23126535356657, 0.01109499526762, 0.01030153766262, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.52895881204315, -0.00000000000000, -0.00000000000000, 0.11271578616360, 0.08945029183212, 0.10777729212931, 0.08553112449288, 0.09384785305309, 0.07447668342215, -0.00310973669075, -0.06323108329847, -0.00302748752528, -0.06155887846385, -0.00313885888624, -0.06382316982296, 0.00015820520724, 0.00300278963678, 0.00586404765501, 0.02939767842756, -0.00016283581422, 0.00536890041324, 0.01109499526762, 0.13938207204499, -0.00030809138916, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.52895881204315, 0.00000000000000, 0.00516592964059, 0.00409963794481, 0.23244284220081, 0.18446462404997, -0.24701084342072, -0.19602524499815, -0.00129823047812, -0.02639725728123, 0.00153171264705, 0.03114480634321, 0.00004393356621, 0.00089331172850, 0.00014689117584, 0.00278804540156, 0.00544468080575, -0.00016283581422, 0.02942186575520, 0.00498494397517, 0.01030153766262, -0.00030809138916, 0.13942783536437, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.52895881204315, 0.27391268082676, 0.21737478014755, -0.14870816683595, -0.11801351174692, -0.13393423370908, -0.10628881960354, 0.00094604275300, 0.01923613285215, 0.00085960816979, 0.01747868964249, -0.00154807357794, -0.03147735098786, 0.03627126005389, 0.21052066164641, 0.17951793448373, -0.00536890041324, -0.00498494397517, 0.43742292602199, 0.56204331679435, -0.01680920964511, -0.01560710050444, 0.18031400188717, 0.76135740608460, 0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.34031690548460, 0.67062731492281, 0.34031695289063, 0.67062737030108, 0.34031737801021, 0.67062786691263, 0.04517041139024, 0.20221562667868, 0.04517032320156, 0.20221539199371, 0.04517043960987, 0.20221570177584, -0.09038250777831, -0.42200122290709, -0.23126535356657, 0.01109499526762, 0.01030153766262, -0.56204331679435, -0.28474473177771, 0.02159809964134, 0.02005351345626, 0.00000000000000, 0.00000000000000, 0.52895881204315, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.19176667863575, 0.22401597325442, 0.18336463168291, 0.21420094658683, 0.15966574552413, 0.18651653837804, -0.11864509679157, -0.31573448993903, -0.11550740200871, -0.30738480224809, -0.11975607686306, -0.31869091315997, 0.00270309863333, 0.01262092585105, 0.01109499526762, 0.13938207204499, -0.00030809138916, 0.01680920964511, 0.02159809964134, 0.43677698473249, -0.00059974685533, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.52895881204315, 0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00878894787465, 0.01026698030142, 0.39546174621222, 0.46196629960800, -0.42024584670042, -0.49091807599058, -0.04953110055667, -0.13181056102010, 0.05843926589535, 0.15551680566247, 0.00167618606718, 0.00446061095493, 0.00250978677370, 0.01171834145538, 0.01030153766262, -0.00030809138916, 0.13942783536437, 0.01560710050444, 0.02005351345626, -0.00059974685533, 0.43686607000996, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.52895881204315, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.46601569155645, 0.54438528861494, -0.25300151545280, -0.29554862134814, -0.22786572713885, -0.26618562736432, 0.03609416010423, 0.09605261016649, 0.03279653693343, 0.08727715145860, -0.05906325359278, -0.15717717810016, 0.00000201298284, 0.00132075867257, 0.00313884064795, -0.00004348617703, 0.00154804262029, 0.04517024141366, 0.11975617197283, -0.00165912790081, 0.05906246256191, 0.01936363367582, 0.19687092023383, 0.11271578616360, 0.00516592964059, 0.27391268082676, 0.34031690548460, 0.19176667863575, 0.00878894787465, 0.46601569155645, 1.00000000000000, 0.64589812339497, 0.00864749381029, 0.12157038108703, 0.00864751028315, 0.12157046545952, 0.00007101679065, 0.02104047963264, 0.00007111844531, 0.02105147921229, 0.00000074574990, 0.00398531987458, 0.00852929429132, 0.06414574506885, 0.06382323888776, -0.00088422095174, 0.03147693847665, 0.20221517434217, 0.31869173496317, -0.00441522420535, 0.15717535351181, 0.08126591223131, 0.40307553361171, 0.08945029183212, 0.00409963794481, 0.21737478014755, 0.67062731492281, 0.22401597325442, 0.01026698030142, 0.54438528861494, 0.64589812339497, 1.00000000000000, 0.12157038108703, 0.36334645092564, 0.12157046545952, 0.36334660395387, 0.02104047963264, 0.12519798346605, 0.02105147921229, 0.12523777123352, 0.00398531987458, 0.04549755009514, 0.00000201299747, 0.00132076396622, 0.00310970417580, 0.00129798720538, -0.00094641848164, 0.04517031385342, 0.11864424188648, 0.04952197998742, -0.03610861256037, 0.01936364597010, 0.19687098662820, 0.10777729212931, 0.23244284220081, -0.14870816683595, 0.34031695289063, 0.18336463168291, 0.39546174621222, -0.25300151545280, 0.00864749381029, 0.12157038108703, 1.00000000000000, 0.64589812339497, 0.00864754110744, 0.12157062333841, 0.00000074575085, 0.00398532171675, 0.00007101641488, 0.02104043895441, 0.00007111967160, 0.02105161184254, 0.00852930957856, 0.06414583914144, 0.06323063664233, 0.02639240027670, -0.01924383791544, 0.20221536711673, 0.31573249216470, 0.13178640538930, -0.09609115496058, 0.08126592485038, 0.40307558443743, 0.08553112449288, 0.18446462404997, -0.11801351174692, 0.67062737030108, 0.21420094658683, 0.46196629960800, -0.29554862134814, 0.12157038108703, 0.36334645092564, 0.64589812339497, 1.00000000000000, 0.12157062333841, 0.36334689030215, 0.00398532171675, 0.04549756289772, 0.02104043895441, 0.12519783630928, 0.02105161184254, 0.12523825093489, 0.00000201305911, 0.00132078627698, 0.00302753686096, -0.00153194938190, -0.00085922913249, 0.04517061915832, 0.11550814510939, -0.05844772157416, -0.03278175225457, 0.01936375622093, 0.19687158202877, 0.09384785305309, -0.24701084342072, -0.13393423370908, 0.34031737801021, 0.15966574552413, -0.42024584670042, -0.22786572713885, 0.00864751028315, 0.12157046545952, 0.00864754110744, 0.12157062333841, 1.00000000000000, 0.64589812339497, 0.00007112075720, 0.02105172925462, 0.00000074576264, 0.00398534470595, 0.00007101868596, 0.02104068480711, 0.00852937400828, 0.06414623562007, 0.06155924802670, -0.03114929934657, -0.01747080273757, 0.20221617958385, 0.30738596067959, -0.15553889319746, -0.08723757445839, 0.08126603801354, 0.40307604022348, 0.07447668342215, -0.19602524499815, -0.10628881960354, 0.67062786691263, 0.18651653837804, -0.49091807599058, -0.26618562736432, 0.12157046545952, 0.36334660395387, 0.12157062333841, 0.36334689030215, 0.64589812339497, 1.00000000000000, 0.02105172925462, 0.12523867559369, 0.00398534470595, 0.04549772266647, 0.02104068480711, 0.12519872569960, 0.01936364917292, 0.19687100392480, -0.10777995984544, -0.23248289019685, 0.14864365644642, 0.34031696524051, -0.18336915941501, -0.39552985745125, 0.25289174698190, 0.00000201301716, 0.00132077109389, -0.00310973669075, -0.00129823047812, 0.00094604275300, 0.04517041139024, -0.11864509679157, -0.04953110055667, 0.03609416010423, 0.00007101679065, 0.02104047963264, 0.00000074575085, 0.00398532171675, 0.00007112075720, 0.02105172925462, 1.00000000000000, 0.64589812339497, 0.00864750286217, 0.12157042744993, 0.00864748844429, 0.12157035360285, 0.08126592813781, 0.40307559767819, -0.08553323626446, -0.18449639434980, 0.11796230954508, 0.67062738472784, -0.21420623256681, -0.46204585821634, 0.29542038879956, 0.00852933016215, 0.06414596580588, -0.06323108329847, -0.02639725728123, 0.01923613285215, 0.20221562667868, -0.31573448993903, -0.13181056102010, 0.09605261016649, 0.02104047963264, 0.12519798346605, 0.00398532171675, 0.04549756289772, 0.02105172925462, 0.12523867559369, 0.64589812339497, 1.00000000000000, 0.12157042744993, 0.36334653501506, 0.12157035360285, 0.36334640107698, 0.01936365384454, 0.19687102915351, -0.09384762989784, 0.24697423962126, 0.13400043121686, 0.34031698325398, -0.15966567009192, 0.42018437234376, 0.22797878503845, 0.00000201299935, 0.00132076464935, -0.00302748752528, 0.00153171264705, 0.00085960816979, 0.04517032320156, -0.11550740200871, 0.05843926589535, 0.03279653693343, 0.00007111844531, 0.02105147921229, 0.00007101641488, 0.02104043895441, 0.00000074576264, 0.00398534470595, 0.00864750286217, 0.12157042744993, 1.00000000000000, 0.64589812339497, 0.00864752726859, 0.12157055245727, 0.08126593293284, 0.40307561699107, -0.07447665397417, 0.19599658515435, 0.10634156407564, 0.67062740577061, -0.18651653871154, 0.49084649633901, 0.26631782436734, 0.00852931155134, 0.06414585128124, -0.06155887846385, 0.03114480634321, 0.01747868964249, 0.20221539199371, -0.30738480224809, 0.15551680566247, 0.08727715145860, 0.02105147921229, 0.12523777123352, 0.02104043895441, 0.12519783630928, 0.00398534470595, 0.04549772266647, 0.12157042744993, 0.36334653501506, 0.64589812339497, 1.00000000000000, 0.12157055245727, 0.36334676174355, 0.01936370483231, 0.19687130450898, -0.11271326806701, -0.00509022391160, -0.27391562534731, 0.34031717985950, -0.19176214056360, -0.00866013602460, -0.46602008398150, 0.00000201302286, 0.00132077315609, -0.00313885888624, 0.00004393356621, -0.00154807357794, 0.04517043960987, -0.11975607686306, 0.00167618606718, -0.05906325359278, 0.00000074574990, 0.00398531987458, 0.00007111967160, 0.02105161184254, 0.00007101868596, 0.02104068480711, 0.00864748844429, 0.12157035360285, 0.00864752726859, 0.12157055245727, 1.00000000000000, 0.64589812339497, 0.08126598526753, 0.40307582777895, -0.08944817024160, -0.00403955295434, -0.21737681737109, 0.67062763543913, -0.22401059818092, -0.01011650290040, -0.54439024027476, 0.00852933611745, 0.06414600245280, -0.06382316982296, 0.00089331172850, -0.03147735098786, 0.20221570177584, -0.31869091315997, 0.00446061095493, -0.15717717810016, 0.00398531987458, 0.04549755009514, 0.02105161184254, 0.12523825093489, 0.02104068480711, 0.12519872569960, 0.12157035360285, 0.36334640107698, 0.12157055245727, 0.36334676174355, 0.64589812339497, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}