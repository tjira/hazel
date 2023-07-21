#include "../include/integral.h"

int test_int_nuclear_ethylene_321g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "3-21G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix V = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Vexp(system.shells.nbf(), system.shells.nbf()); Vexp << -38.35278658340739, -3.05860490842776, -0.06771851259754, -0.00094034169179, -0.01792900000415, -4.17944916512248, -0.01004714822188, -0.00013951506013, -0.00266006020514, -0.00000026205916, -0.17181899926011, 0.46820361453741, 0.00649568403910, 0.12401876080673, -1.30189003467416, 2.69046035306431, 0.03732644871214, 0.71265482189246, -0.42071389271324, -1.87621750164862, -0.42071058232354, -1.87621399338494, -0.00008036913491, -0.24015017873438, -0.00008035975124, -0.24014337921706, -3.05860490842776, -10.15753787235205, -0.49807421374311, -0.00691625128057, -0.13186932969837, -7.27770728921757, -0.31841081458391, -0.00442131606922, -0.08430347485697, -0.17181900013956, -1.32516604885449, 2.08812696507705, 0.02896966010783, 0.55311084467654, -2.95783835479611, 4.83819859852170, 0.06712254491396, 1.28156294485716, -1.69025628641932, -3.69835882871720, -1.69024615187770, -3.69834979159382, -0.02060593640279, -0.73951613139317, -0.02060476047446, -0.73949935271563, -0.06771851259754, -0.49807421374311, -10.87099837991164, 0.00183537985082, -0.09605784211746, -0.49003772054998, -4.72092840847968, 0.00117814403478, -0.08994223544426, -0.46820364852627, -2.08812876216158, 2.83115735395963, 0.04915890444643, 0.93454169518005, -2.19367665994200, 1.77450458889720, 0.04963804319416, 0.92130394614721, 1.11356710766407, 0.70753974061583, 1.32570821538409, 0.87566138240934, -0.04762239223515, -0.72774474851187, -0.04555930821755, -0.69683238291033, -0.00094034169179, -0.00691625128057, 0.00183537985082, -10.63756675374400, -0.02047409411844, -0.00680447367908, 0.00117814403478, -4.49219425922328, -0.01766426802037, -0.00649569604663, -0.02897021534458, 0.04915974303481, -0.70026277319498, 0.01237584943379, -0.03043568442122, 0.04964099770177, -1.72900235293430, 0.00892362215461, 2.05963059957389, 1.62976974279835, -2.02587810231513, -1.60788608254186, 0.01920221079730, 0.28774150015103, -0.02049476325466, -0.30751356512790, -0.01792900000415, -0.13186932969837, -0.09605784211746, -0.02047409411844, -10.53276192210044, -0.12974373224643, -0.08994223544426, -0.01766426802037, -4.40430620826497, -0.12401864272687, -0.55310401186265, 0.93453024681680, 0.01237546847635, -0.44938930636463, -0.58104786284875, 0.92126691447298, 0.00892232579460, -1.45939961880275, 0.61688423052874, 0.44254815429929, 0.03015463076671, -0.02241783163117, -0.00948628234543, -0.14585987752458, -0.01518679083348, -0.23134502639893, -4.17944916512248, -7.27770728921757, -0.49003772054998, -0.00680447367908, -0.12974373224643, -8.01295300341450, -0.81631550350901, -0.01133332251379, -0.21614807677428, -1.30189004779992, -2.95783835671392, 2.19366797081048, 0.03043259371206, 0.58108060992910, -4.63719772678834, 4.89544812010701, 0.06791485646344, 1.29675171332459, -2.41854269851811, -4.74219095792144, -2.41853241502883, -4.74217728009383, -0.41603577422746, -1.77082096473008, -0.41602587851398, -1.77079548592316, -0.01004714822188, -0.31841081458391, -4.72092840847968, 0.00117814403478, -0.08994223544426, -0.81631550350901, -6.91858150274687, -0.00107475492987, -0.25393508034358, -2.69046056740117, -4.83820422554467, 1.77448993060462, 0.04963467247104, 0.92132947548430, -4.89546578873218, 0.89702864083047, 0.06472127282185, 1.14697638242997, 1.50882545882410, 1.23407500141410, 1.79484481887313, 1.52407940385226, -1.03514112970226, -2.60368364103879, -0.98734503705537, -2.49913292650047, -0.00013951506013, -0.00442131606922, 0.00117814403478, -4.49219425922328, -0.01766426802037, -0.01133332251379, -0.00107475492987, -6.19001358579711, -0.03760631872779, -0.03732652766261, -0.06712451350488, 0.04963798911582, -1.72900235220209, 0.00892324628083, -0.06792077999570, 0.06472726747300, -3.51963207780390, 0.00293869412413, 2.77696356268248, 2.81147238863533, -2.73125558057339, -2.77334434459356, 0.44615081511358, 0.97139948593612, -0.47422131669848, -1.04222974614493, -0.00266006020514, -0.08430347485697, -0.08994223544426, -0.01766426802037, -4.40430620826497, -0.21614807677428, -0.25393508034358, -0.03760631872779, -6.02529345950775, -0.71265412908348, -1.28154174552875, 0.92129250588095, 0.00892185485837, -1.45938506712197, -1.29668488065963, 1.14690088690776, 0.00293605885124, -3.12842282272704, 0.83368272292503, 0.76697043781026, 0.04263525819456, -0.03507536863936, -0.20166837996008, -0.53099502792617, -0.33384390914026, -0.82017940938792, -0.00000026205916, -0.17181900013956, -0.46820364852627, -0.00649569604663, -0.12401864272687, -1.30189004779992, -2.69046056740117, -0.03732652766261, -0.71265412908348, -38.35278687821926, -3.05860496486874, 0.06771430287084, 0.00093849711752, 0.01794469283227, -4.17944921828119, 0.01004652364012, 0.00013924138738, 0.00266238849269, -0.00008037253772, -0.24015264617550, -0.00008035921084, -0.24014299050314, -0.42071028981424, -1.87621370742943, -0.42071367245666, -1.87621728034415, -0.17181899926011, -1.32516604885449, -2.08812876216158, -0.02897021534458, -0.55310401186265, -2.95783835671392, -4.83820422554467, -0.06712451350488, -1.28154174552875, -3.05860496486874, -10.15753816784397, 0.49804341460642, 0.00690275947288, 0.13198414037097, -7.27770751733226, 0.31839189336994, 0.00441304386873, 0.08437401075392, -0.02060638687143, -0.73952220729784, -0.02060467409492, -0.73949834213934, -1.69024766224528, -3.69835207712917, -1.69025337549739, -3.69835548241257, 0.46820361453741, 2.08812696507705, 2.83115735395963, 0.04915974303481, 0.93453024681680, 2.19366797081048, 1.77448993060462, 0.04963798911582, 0.92129250588095, 0.06771430287084, 0.49804341460642, -10.87101393242901, 0.00182871079100, -0.09603156783323, 0.49000848275285, -4.72094130948083, 0.00117249938734, -0.08992031003302, 0.04556016300401, 0.69681176171817, 0.04761771524920, 0.72769933922102, -1.32589485076723, -0.87581833260505, -1.11384185920320, -0.70777023954669, 0.00649568403910, 0.02896966010783, 0.04915890444643, -0.70026277319498, 0.01237546847635, 0.03043259371206, 0.04963467247104, -1.72900235220209, 0.00892185485837, 0.00093849711752, 0.00690275947288, 0.00182871079100, -10.63756679209036, -0.02047539187690, 0.00679168853954, 0.00117249938734, -4.49219420700376, -0.01766537296746, 0.02049546986081, 0.30750822744192, -0.01920201203332, -0.28775023284578, 2.02577195276501, 1.60779693801960, -2.05972957230264, -1.62985408489754, 0.12401876080673, 0.55311084467654, 0.93454169518005, 0.01237584943379, -0.44938930636463, 0.58108060992910, 0.92132947548430, 0.00892324628083, -1.45938506712197, 0.01794469283227, 0.13198414037097, -0.09603156783323, -0.02047539187690, -10.53274722078739, 0.12985272674454, -0.08992031003302, -0.01766537296746, -4.40429384867396, 0.01519616976068, 0.23147989200069, 0.00949413471440, 0.14598779398723, -0.02928731172943, 0.02314439576201, -0.61603960905318, -0.44184069531858, -1.30189003467416, -2.95783835479611, -2.19367665994200, -0.03043568442122, -0.58104786284875, -4.63719772678834, -4.89546578873218, -0.06792077999570, -1.29668488065963, -4.17944921828119, -7.27770751733226, 0.49000848275285, 0.00679168853954, 0.12985272674454, -8.01295332869018, 0.81627590439882, 0.01131629685632, 0.21629588532227, -0.41603977480708, -1.77083048339339, -0.41602495162569, -1.77079359184481, -2.41853603286509, -4.74218387701585, -2.41853778218763, -4.74218283048825, 2.69046035306431, 4.83819859852170, 1.77450458889720, 0.04964099770177, 0.92126691447298, 4.89544812010701, 0.89702864083047, 0.06472726747300, 1.14690088690776, 0.01004652364012, 0.31839189336994, -4.72094130948083, 0.00117249938734, -0.08992031003302, 0.81627590439882, -6.91860101912823, -0.00108485779949, -0.25390232862522, 0.98731594970681, 2.49903799113040, 1.03507235698558, 2.60355132452540, -1.79509715043664, -1.52434543871190, -1.50919468447601, -1.23446957783122, 0.03732644871214, 0.06712254491396, 0.04963804319416, -1.72900235293430, 0.00892232579460, 0.06791485646344, 0.06472127282185, -3.51963207780390, 0.00293605885124, 0.00013924138738, 0.00441304386873, 0.00117249938734, -4.49219420700376, -0.01766537296746, 0.01131629685632, -0.00108485779949, -6.19001354590425, -0.03760841067756, 0.47421480222389, 1.04219987304603, -0.44616288125620, -0.97143761727271, 2.73111485213277, 2.77319471880477, -2.77709561356151, -2.81161373363532, 0.71265482189246, 1.28156294485716, 0.92130394614721, 0.00892362215461, -1.45939961880275, 1.29675171332459, 1.14697638242997, 0.00293869412413, -3.12842282272704, 0.00266238849269, 0.08437401075392, -0.08992031003302, -0.01766537296746, -4.40429384867396, 0.21629588532227, -0.25390232862522, -0.03760841067756, -6.02527501627952, 0.33404480031096, 0.82062494715538, 0.20185812523004, 0.53142977683217, -0.04146828652689, 0.03631358019498, -0.83254601218203, -0.76576562001294, -0.42071389271324, -1.69025628641932, 1.11356710766407, 2.05963059957389, 0.61688423052874, -2.41854269851811, 1.50882545882410, 2.77696356268248, 0.83368272292503, -0.00008037253772, -0.02060638687143, 0.04556016300401, 0.02049546986081, 0.01519616976068, -0.41603977480708, 0.98731594970681, 0.47421480222389, 0.33404480031096, -6.84873559982649, -4.08837166002777, -0.05001909390064, -0.75646552943263, -0.00078208439376, -0.17220639778775, -0.00000870894514, -0.03440920630803, -1.87621750164862, -3.69835882871720, 0.70753974061583, 1.62976974279835, 0.44254815429929, -4.74219095792144, 1.23407500141410, 2.81147238863533, 0.76697043781026, -0.24015264617550, -0.73952220729784, 0.69681176171817, 0.30750822744192, 0.23147989200069, -1.77083048339339, 2.49903799113040, 1.04219987304603, 0.82062494715538, -4.08837166002777, -5.57432435143988, -0.75646392240359, -2.22199068151056, -0.17220631595665, -0.90833879778847, -0.03440928670998, -0.38187024952577, -0.42071058232354, -1.69024615187770, 1.32570821538409, -2.02587810231513, 0.03015463076671, -2.41853241502883, 1.79484481887313, -2.73125558057339, 0.04263525819456, -0.00008035921084, -0.02060467409492, 0.04761771524920, -0.01920201203332, 0.00949413471440, -0.41602495162569, 1.03507235698558, -0.44616288125620, 0.20185812523004, -0.05001909390064, -0.75646392240359, -6.84871740984444, -4.08835998736050, -0.00000870863213, -0.03440880859771, -0.00078185687959, -0.17218817393779, -1.87621399338494, -3.69834979159382, 0.87566138240934, -1.60788608254186, -0.02241783163117, -4.74217728009383, 1.52407940385226, -2.77334434459356, -0.03507536863936, -0.24014299050314, -0.73949834213934, 0.72769933922102, -0.28775023284578, 0.14598779398723, -1.77079359184481, 2.60355132452540, -0.97143761727271, 0.53142977683217, -0.75646552943263, -2.22199068151056, -4.08835998736050, -5.57430760772912, -0.03440874703409, -0.38186720028663, -0.17218827522691, -0.90828259223592, -0.00008036913491, -0.02060593640279, -0.04762239223515, 0.01920221079730, -0.00948628234543, -0.41603577422746, -1.03514112970226, 0.44615081511358, -0.20166837996008, -0.42071028981424, -1.69024766224528, -1.32589485076723, 2.02577195276501, -0.02928731172943, -2.41853603286509, -1.79509715043664, 2.73111485213277, -0.04146828652689, -0.00078208439376, -0.17220631595665, -0.00000870863213, -0.03440874703409, -6.84872911785373, -4.08836755863680, -0.05001977167716, -0.75646872153794, -0.24015017873438, -0.73951613139317, -0.72774474851187, 0.28774150015103, -0.14585987752458, -1.77082096473008, -2.60368364103879, 0.97139948593612, -0.53099502792617, -1.87621370742943, -3.69835207712917, -0.87581833260505, 1.60779693801960, 0.02314439576201, -4.74218387701585, -1.52434543871190, 2.77319471880477, 0.03631358019498, -0.17220639778775, -0.90833879778847, -0.03440880859771, -0.38186720028663, -4.08836755863680, -5.57431938002887, -0.75646833060329, -2.22199658342351, -0.00008035975124, -0.02060476047446, -0.04555930821755, -0.02049476325466, -0.01518679083348, -0.41602587851398, -0.98734503705537, -0.47422131669848, -0.33384390914026, -0.42071367245666, -1.69025337549739, -1.11384185920320, -2.05972957230264, -0.61603960905318, -2.41853778218763, -1.50919468447601, -2.77709561356151, -0.83254601218203, -0.00000870894514, -0.03440928670998, -0.00078185687959, -0.17218827522691, -0.05001977167716, -0.75646833060329, -6.84872190876833, -4.08836282006906, -0.24014337921706, -0.73949935271563, -0.69683238291033, -0.30751356512790, -0.23134502639893, -1.77079548592316, -2.49913292650047, -1.04222974614493, -0.82017940938792, -1.87621728034415, -3.69835548241257, -0.70777023954669, -1.62985408489754, -0.44184069531858, -4.74218283048825, -1.23446957783122, -2.81161373363532, -0.76576562001294, -0.03440920630803, -0.38187024952577, -0.17218817393779, -0.90828259223592, -0.75646872153794, -2.22199658342351, -4.08836282006906, -5.57431080259266;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR: " << V << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR NORM: " << V.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR: " << Vexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR NORM: " << Vexp.norm() << std::endl;

    // return success or failure based on the error
    return (V - Vexp).norm() > 1e-8;
}
