#include "../include/integral.h"

int test_integral_overlap_ammonia_ccpvdz(int, char**) {
    // initialize the system
    System system("../example/molecule/ammonia.xyz", "CC-PVDZ", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, -0.00000016705639, 0.18957710152340, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.06964491584187, 0.07822966134016, -0.03382759201800, -0.07395970220641, -0.06229946440363, 0.06964494092516, 0.07822967471323, -0.08495302264535, 0.03623667731133, 0.04433361173749, 0.06964493735224, 0.07822967280834, 0.04791444667254, -0.05050212309729, 0.07516196624946, -0.00000016705639, 1.00000000000000, 0.93585039893372, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.53927357250987, 0.58498469641444, -0.12416930260839, -0.27148029452423, -0.22867962472048, 0.53927367926242, 0.58498476886552, -0.31183276822150, 0.13301214065466, 0.16273314877886, 0.53927366405637, 0.58498475854544, 0.17587715271182, -0.18537560658806, 0.27589325421037, 0.18957710152340, 0.93585039893372, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.62624709722769, 0.69069837430388, -0.10566673920995, -0.23102680678335, -0.19460389774578, 0.62624719862714, 0.69069845246860, -0.26536629963979, 0.11319188735041, 0.13848414253086, 0.62624718418359, 0.69069844133466, 0.14966954937255, -0.15775263059983, 0.23478216696090, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.81606765954104, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.14756134772906, 0.10175302154383, 0.22596119121409, -0.14586482875896, -0.12286827063811, 0.37057838473002, 0.25553754896432, -0.12808999913427, 0.17947789498680, 0.21958148213880, -0.20901033695648, -0.14412602697028, 0.15882724011248, 0.14107821614911, -0.20996574936672, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.81606765954104, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.32262400851377, 0.22246996380878, -0.14586482875896, -0.02623813070980, -0.26863575456905, -0.15807005952072, -0.10899943774108, 0.17947789498680, 0.21612054789971, -0.09366239203783, 0.22029818767857, 0.15190972370653, 0.14107821614911, 0.14397944248877, 0.22130519826727, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.81606765954104, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.27176019284214, 0.18739609784396, -0.12286827063811, -0.26863575456905, 0.06639311075661, -0.19339015511566, -0.13335490754024, 0.21958148213880, -0.09366239203783, 0.17808590818806, -0.32786829407578, -0.22608620837976, -0.20996574936672, 0.22130519826727, -0.03669026634241, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81606765954104, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.19329054908506, 0.15043565212239, 0.29623432124107, -0.10497410103949, -0.08842423746869, 0.48542043292638, 0.37779671540357, 0.04143544011116, 0.12916428493801, 0.15802550574917, -0.27378253475146, -0.21308156047834, 0.24792025248385, 0.10152931375631, -0.15110538697911, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81606765954104, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.42260505690251, 0.32890830734499, -0.10497410103949, 0.11473479646203, -0.19332828264962, -0.20705588854339, -0.16114903553828, 0.12916428493801, 0.28915234797346, -0.06740571530572, 0.28856848470777, 0.22458928244719, 0.10152931375631, 0.23723478491361, 0.15926601231641, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81606765954104, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.35597856554122, 0.27705373026397, -0.08842423746869, -0.19332828264962, 0.18139844241315, -0.25332166334616, -0.19715711548397, 0.15802550574917, -0.06740571530572, 0.26178007997583, -0.42947451271462, -0.33425470122841, -0.15110538697911, 0.15926601231641, 0.10721278679878, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.08483142264982, 0.01842903664506, 0.18786397578858, -0.06647611224486, -0.16233343564030, -0.10437996269342, -0.02267581342520, 0.13714132907588, 0.20085378441554, -0.14214015991307, -0.08204764664483, -0.01782427554227, 0.06771420885650, -0.05153589834552, -0.18942219414421, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.15623199525609, 0.03394031451754, -0.16233343564030, -0.12242746058561, -0.02295721168857, 0.05447176092300, 0.01183360729166, -0.14214015991307, -0.10481762057435, -0.06105340029942, -0.12870570101982, -0.02796041047672, -0.18942219414421, -0.08084282969089, -0.10867339330240, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.01123778555168, 0.00244133076250, -0.08456156457582, -0.18488304254134, 0.24695611890199, -0.04501727526387, -0.00977968671815, -0.06557028967551, 0.02796898043406, -0.25234485355821, 0.06315672159533, 0.01372035462436, 0.19618709343382, -0.20678240970635, -0.17807872671232, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.07145718598126, 0.01552357673450, 0.15824597345866, -0.16233343564030, -0.01050013950441, -0.12770323005794, -0.02774262936934, 0.16778498713271, -0.14214015991307, 0.14313318115926, 0.12211095434713, 0.02652774803443, -0.10077859639437, -0.18942219414421, 0.10310508130679, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.07333648368887, -0.01593184109267, 0.20244089895330, -0.10940563566237, 0.14033671702214, 0.10009222601126, 0.02174433276059, 0.05585046992871, 0.24663842452299, 0.13630130385353, -0.00431755832380, -0.00093795925149, -0.18516503396432, -0.18177003756403, -0.00996788335174, 0.06964491584187, 0.53927357250987, 0.62624709722769, 0.14756134772906, 0.32262400851377, 0.27176019284214, 0.19329054908506, 0.42260505690251, 0.35597856554122, 0.08483142264982, 0.15623199525609, 0.01123778555168, 0.07145718598126, -0.07333648368887, 1.00000000000000, 0.90368690205096, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.33143465720527, 0.44049315146007, -0.04812691684372, 0.10373347208569, 0.10037914494187, 0.33143462787800, 0.44049312096743, 0.07694792422622, 0.02208185543160, 0.12939941607376, 0.07822966134016, 0.58498469641444, 0.69069837430388, 0.10175302154383, 0.22246996380878, 0.18739609784396, 0.15043565212239, 0.32890830734499, 0.27705373026397, 0.01842903664506, 0.03394031451754, 0.00244133076250, 0.01552357673450, -0.01593184109267, 0.90368690205096, 1.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.44049315146007, 0.56094958715653, -0.05209446268376, 0.11228517937629, 0.10865432409440, 0.44049312096743, 0.56094955525101, 0.08329145812238, 0.02390226839047, 0.14006701484606, -0.03382759201800, -0.12416930260839, -0.10566673920995, 0.22596119121409, -0.14586482875896, -0.12286827063811, 0.29623432124107, -0.10497410103949, -0.08842423746869, 0.18786397578858, -0.16233343564030, -0.08456156457582, 0.15824597345866, 0.20244089895330, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.04812691684372, 0.05209446268376, 0.00991280613426, 0.04739822502907, 0.04586555529781, -0.07694792422622, -0.08329145812238, -0.02431147898825, -0.01613197719864, -0.09453319881048, -0.07395970220641, -0.27148029452423, -0.23102680678335, -0.14586482875896, -0.02623813070980, -0.26863575456905, -0.10497410103949, 0.11473479646203, -0.19332828264962, -0.06647611224486, -0.12242746058561, -0.18488304254134, -0.16233343564030, -0.10940563566237, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.10373347208569, -0.11228517937629, 0.04739822502907, -0.07025973024210, -0.09885929978914, -0.02208185543160, -0.02390226839047, -0.01613197719864, 0.02727368138194, -0.02712832673019, -0.06229946440363, -0.22867962472048, -0.19460389774578, -0.12286827063811, -0.26863575456905, 0.06639311075661, -0.08842423746869, -0.19332828264962, 0.18139844241315, -0.16233343564030, -0.02295721168857, 0.24695611890199, -0.01050013950441, 0.14033671702214, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.10037914494187, -0.10865432409440, 0.04586555529781, -0.09885929978914, -0.06375947528286, -0.12939941607376, -0.14006701484606, -0.09453319881048, -0.02712832673019, -0.12706858134654, 0.06964494092516, 0.53927367926242, 0.62624719862714, 0.37057838473002, -0.15807005952072, -0.19339015511566, 0.48542043292638, -0.20705588854339, -0.25332166334616, -0.10437996269342, 0.05447176092300, -0.04501727526387, -0.12770323005794, 0.10009222601126, 0.33143465720527, 0.44049315146007, 0.04812691684372, -0.10373347208569, -0.10037914494187, 1.00000000000000, 0.90368690205096, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33143367268533, 0.44049212781785, 0.12507411910976, -0.08165113267623, 0.02902012315809, 0.07822967471323, 0.58498476886552, 0.69069845246860, 0.25553754896432, -0.10899943774108, -0.13335490754024, 0.37779671540357, -0.16114903553828, -0.19715711548397, -0.02267581342520, 0.01183360729166, -0.00977968671815, -0.02774262936934, 0.02174433276059, 0.44049315146007, 0.56094958715653, 0.05209446268376, -0.11228517937629, -0.10865432409440, 0.90368690205096, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.44049212781785, 0.56094851608334, 0.13538548620342, -0.08838261964276, 0.03141260167495, -0.08495302264535, -0.31183276822150, -0.26536629963979, -0.12808999913427, 0.17947789498680, 0.21958148213880, 0.04143544011116, 0.12916428493801, 0.15802550574917, 0.13714132907588, -0.14214015991307, -0.06557028967551, 0.16778498713271, 0.05585046992871, -0.04812691684372, -0.05209446268376, 0.00991280613426, 0.04739822502907, 0.04586555529781, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.12507411910976, -0.13538548620342, -0.11661905131911, 0.09695829191675, -0.03446053325162, 0.03623667731133, 0.13301214065466, 0.11319188735041, 0.17947789498680, 0.21612054789971, -0.09366239203783, 0.12916428493801, 0.28915234797346, -0.06740571530572, 0.20085378441554, -0.10481762057435, 0.02796898043406, -0.14214015991307, 0.24663842452299, 0.10373347208569, 0.11228517937629, 0.04739822502907, -0.07025973024210, -0.09885929978914, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.08165113267623, 0.08838261964276, 0.09695829191675, -0.03139375733780, 0.02249659316131, 0.04433361173749, 0.16273314877886, 0.13848414253086, 0.21958148213880, -0.09366239203783, 0.17808590818806, 0.15802550574917, -0.06740571530572, 0.26178007997583, -0.14214015991307, -0.06105340029942, -0.25234485355821, 0.14313318115926, 0.13630130385353, 0.10037914494187, 0.10865432409440, 0.04586555529781, -0.09885929978914, -0.06375947528286, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.02902012315809, -0.03141260167495, -0.03446053325162, 0.02249659316131, 0.02390709536883, 0.06964493735224, 0.53927366405637, 0.62624718418359, -0.20901033695648, 0.22029818767857, -0.32786829407578, -0.27378253475146, 0.28856848470777, -0.42947451271462, -0.08204764664483, -0.12870570101982, 0.06315672159533, 0.12211095434713, -0.00431755832380, 0.33143462787800, 0.44049312096743, -0.07694792422622, -0.02208185543160, -0.12939941607376, 0.33143367268533, 0.44049212781785, -0.12507411910976, 0.08165113267623, -0.02902012315809, 1.00000000000000, 0.90368690205096, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.07822967280834, 0.58498475854544, 0.69069844133466, -0.14412602697028, 0.15190972370653, -0.22608620837976, -0.21308156047834, 0.22458928244719, -0.33425470122841, -0.01782427554227, -0.02796041047672, 0.01372035462436, 0.02652774803443, -0.00093795925149, 0.44049312096743, 0.56094955525101, -0.08329145812238, -0.02390226839047, -0.14006701484606, 0.44049212781785, 0.56094851608334, -0.13538548620342, 0.08838261964276, -0.03141260167495, 0.90368690205096, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.04791444667254, 0.17587715271182, 0.14966954937255, 0.15882724011248, 0.14107821614911, -0.20996574936672, 0.24792025248385, 0.10152931375631, -0.15110538697911, 0.06771420885650, -0.18942219414421, 0.19618709343382, -0.10077859639437, -0.18516503396432, 0.07694792422622, 0.08329145812238, -0.02431147898825, -0.01613197719864, -0.09453319881048, 0.12507411910976, 0.13538548620342, -0.11661905131911, 0.09695829191675, -0.03446053325162, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.05050212309729, -0.18537560658806, -0.15775263059983, 0.14107821614911, 0.14397944248877, 0.22130519826727, 0.10152931375631, 0.23723478491361, 0.15926601231641, -0.05153589834552, -0.08084282969089, -0.20678240970635, -0.18942219414421, -0.18177003756403, 0.02208185543160, 0.02390226839047, -0.01613197719864, 0.02727368138194, -0.02712832673019, -0.08165113267623, -0.08838261964276, 0.09695829191675, -0.03139375733780, 0.02249659316131, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.07516196624946, 0.27589325421037, 0.23478216696090, -0.20996574936672, 0.22130519826727, -0.03669026634241, -0.15110538697911, 0.15926601231641, 0.10721278679878, -0.18942219414421, -0.10867339330240, -0.17807872671232, 0.10310508130679, -0.00996788335174, 0.12939941607376, 0.14006701484606, -0.09453319881048, -0.02712832673019, -0.12706858134654, 0.02902012315809, 0.03141260167495, -0.03446053325162, 0.02249659316131, 0.02390709536883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}