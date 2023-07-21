#include "../include/integral.h"

int test_int_nuclear_ammonia_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/ammonia.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix V = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Vexp(system.shells.nbf(), system.shells.nbf()); Vexp << -48.36181331889841, -3.68281816961500, -0.01494281122514, -0.01860312213203, 0.01206034377081, -4.42448212422103, -0.00218934794722, -0.00272563887999, 0.00176702284995, -0.77649733027687, 0.00048173788288, -0.00031231086280, -0.77628454008842, -0.00038880914576, -0.77663220589919, -0.77877396359354, -1.98052929171090, -0.77877463611774, -1.98052973866240, -0.77877453936524, -1.98052967432655, -3.68281816961500, -9.47311360487115, -0.10363062480575, -0.12901542793456, 0.08364028150688, -6.59904556032478, -0.06651527211222, -0.08280850001492, 0.05368447436469, -5.02904329888235, 0.02323842939218, -0.01506548516353, -5.01877856965437, -0.01875566700918, -5.03554953117525, -1.96906240428549, -3.41387415861773, -1.96906342194712, -3.41387471002093, -1.96906325284787, -3.41387460963291, -0.01494281122514, -0.10363062480575, -10.27384104458019, 0.02598392166353, -0.01684538866488, -0.09728111401597, -3.76503886816826, 0.02288827521748, -0.01483848447895, -0.15313905984829, -0.04322889709667, 0.06078622452178, -0.05060414229045, -0.03371993981147, -0.04016539104034, -0.78038817614674, -0.49986418239229, -1.94070956584960, -1.18516415332218, 1.07479560916806, 0.59583051202318, -0.01860312213203, -0.12901542793456, 0.02598392166353, -10.26236359379426, -0.02097154492871, -0.12111057839866, 0.02288827521748, -3.75492881055014, -0.01847306503412, -0.04322889709667, -0.05060414229045, -0.03371993981147, -0.17722986427199, 0.02437327903379, -0.08319637930862, -1.69430911602622, -1.04918010353303, 0.80666282656641, 0.42792518913966, -1.16192515382086, -0.73474724178591, 0.01206034377081, 0.08364028150688, -0.01684538866488, -0.02097154492871, -10.28111595187669, 0.07851558579626, -0.01483848447895, -0.01847306503412, -3.77144706946433, 0.06078622452178, -0.03371993981147, -0.04016539104034, 0.02437327903379, -0.08319637930862, 0.11169912557084, -1.40371273934156, -0.79760810057799, 1.01638779502608, 0.63173340668723, 1.71605543999335, 1.04496582891793, -4.42448212422103, -6.59904556032478, -0.09728111401597, -0.12111057839866, 0.07851558579626, -6.56990254944830, -0.12672855201733, -0.15777132145105, 0.10228258600702, -4.56096690114715, 0.03156572586981, -0.02046408170102, -4.54702389741668, -0.02547661203194, -4.56980459311544, -2.49402999229564, -4.07384191008223, -2.49403060550027, -4.07384225555439, -2.49403047914553, -4.07384215899660, -0.00218934794722, -0.06651527211222, -3.76503886816826, 0.02288827521748, -0.01483848447895, -0.12672855201733, -4.71159774734070, 0.04557511388041, -0.02954639859633, -0.15433964189862, -0.03749357710941, 0.06711737170994, -0.05086869359625, -0.04406341492991, -0.03722789098307, -1.00637433263222, -0.84784449325277, -2.50727954094841, -2.02434966446395, 1.39335408056559, 1.03321589060410, -0.00272563887999, -0.08280850001492, 0.02288827521748, -3.75492881055014, -0.01847306503412, -0.15777132145105, 0.04557511388041, -4.69146662722152, -0.03678358165489, -0.03749357710941, -0.05086869359625, -0.04406341492991, -0.17460760512071, 0.01953487918023, -0.08972094677369, -2.18780099911623, -1.78837057562348, 1.04727138200793, 0.74748569588286, -1.49914790900268, -1.24856025157843, 0.00176702284995, 0.05368447436469, -0.01483848447895, -0.01847306503412, -3.77144706946433, 0.10228258600702, -0.02954639859633, -0.03678358165489, -4.72435777227957, 0.06711737170994, -0.04406341492991, -0.03722789098307, 0.01953487918023, -0.08972094677369, 0.10901800307903, -1.81822366463320, -1.37763382405749, 1.31223932619573, 1.07622271835532, 2.21727754036226, 1.78564946328819, -0.77649733027687, -5.02904329888235, -0.15313905984829, -0.04322889709667, 0.06078622452178, -4.56096690114715, -0.15433964189862, -0.03749357710941, 0.06711737170994, -6.84880301407903, 0.04325249501292, 0.00050256845999, -2.28155226264167, -0.01290464713236, -2.29546466453049, -1.19638628901175, -2.50164731305136, -2.58591745914801, -2.79681016419986, -1.44115527671286, -2.53589733618375, 0.00048173788288, 0.02323842939218, -0.04322889709667, -0.05060414229045, -0.03371993981147, 0.03156572586981, -0.03749357710941, -0.05086869359625, -0.04406341492991, 0.04325249501292, -2.28155226264167, -0.01290464713236, 0.01844177374116, -0.01842504427767, 0.01590990373004, -0.57779767994930, -0.11575847404712, 0.69421076853699, 0.14983980931762, 0.55268896823380, 0.12848760608770, -0.00031231086280, -0.01506548516353, 0.06078622452178, -0.03371993981147, -0.04016539104034, -0.02046408170102, 0.06711737170994, -0.04406341492991, -0.03722789098307, 0.00050256845999, -0.01290464713236, -2.29546466453049, -0.01842504427767, 0.01590990373004, -0.03238835373006, -0.48361671303217, -0.11185654757626, 0.86539450582276, 0.17945287932811, -0.81555731740831, -0.17299011941752, -0.77628454008842, -5.01877856965437, -0.05060414229045, -0.17722986427199, 0.02437327903379, -4.54702389741668, -0.05086869359625, -0.17460760512071, 0.01953487918023, -2.28155226264167, 0.01844177374116, -0.01842504427767, -6.81189573159573, -0.03274805180237, -2.29809313771705, -2.18966121336130, -2.71062315079739, -1.21757649135085, -2.48430322034137, -1.52067005632522, -2.56761966339136, -0.00038880914576, -0.01875566700918, -0.03371993981147, 0.02437327903379, -0.08319637930862, -0.02547661203194, -0.04406341492991, 0.01953487918023, -0.08972094677369, -0.01290464713236, -0.01842504427767, 0.01590990373004, -0.03274805180237, -2.29809313771705, -0.01698140645916, -1.04999999203176, -0.22866753733948, -0.36397677357857, -0.08359516632198, 0.87394619642929, 0.18105360357413, -0.77663220589919, -5.03554953117525, -0.04016539104034, -0.08319637930862, 0.11169912557084, -4.56980459311544, -0.03722789098307, -0.08972094677369, 0.10901800307903, -2.29546466453049, 0.01590990373004, -0.03238835373006, -2.29809313771705, -0.01698140645916, -6.85398955043726, -1.80133766769371, -2.61331928146568, -1.38389304656510, -2.54447727121824, -2.22556132251642, -2.72207344438070, -0.77877396359354, -1.96906240428549, -0.78038817614674, -1.69430911602622, -1.40371273934156, -2.49402999229564, -1.00637433263222, -2.18780099911623, -1.81822366463320, -1.19638628901175, -0.57779767994930, -0.48361671303217, -2.18966121336130, -1.04999999203176, -1.80133766769371, -5.88618578836483, -3.50329718115072, -0.24168211108911, -1.15579278567694, -0.24168204120937, -1.15579265647387, -1.98052929171090, -3.41387415861773, -0.49986418239229, -1.04918010353303, -0.79760810057799, -4.07384191008223, -0.84784449325277, -1.78837057562348, -1.37763382405749, -2.50164731305136, -0.11575847404712, -0.11185654757626, -2.71062315079739, -0.22866753733948, -2.61331928146568, -3.50329718115072, -4.44934376284540, -1.15579282491774, -2.44968770116712, -1.15579267660536, -2.44968753663272, -0.77877463611774, -1.96906342194712, -1.94070956584960, 0.80666282656641, 1.01638779502608, -2.49403060550027, -2.50727954094841, 1.04727138200793, 1.31223932619573, -2.58591745914801, 0.69421076853699, 0.86539450582276, -1.21757649135085, -0.36397677357857, -1.38389304656510, -0.24168211108911, -1.15579282491774, -5.88618593288206, -3.50329724085147, -0.24167995581118, -1.15578869035977, -1.98052973866240, -3.41387471002093, -1.18516415332218, 0.42792518913966, 0.63173340668723, -4.07384225555439, -2.02434966446395, 0.74748569588286, 1.07622271835532, -2.79681016419986, 0.14983980931762, 0.17945287932811, -2.48430322034137, -0.08359516632198, -2.54447727121824, -1.15579278567694, -2.44968770116712, -3.50329724085147, -4.44934364606709, -1.15578867125050, -2.44968286903178, -0.77877453936524, -1.96906325284787, 1.07479560916806, -1.16192515382086, 1.71605543999335, -2.49403047914553, 1.39335408056559, -1.49914790900268, 2.21727754036226, -1.44115527671286, 0.55268896823380, -0.81555731740831, -1.52067005632522, 0.87394619642929, -2.22556132251642, -0.24168204120937, -1.15579267660536, -0.24167995581118, -1.15578867125050, -5.88618582218774, -3.50329717311614, -1.98052967432655, -3.41387460963291, 0.59583051202318, -0.73474724178591, 1.04496582891793, -4.07384215899660, 1.03321589060410, -1.24856025157843, 1.78564946328819, -2.53589733618375, 0.12848760608770, -0.17299011941752, -2.56761966339136, 0.18105360357413, -2.72207344438070, -1.15579265647387, -2.44968753663272, -1.15578869035977, -2.44968286903178, -3.50329717311614, -4.44934358216004;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR: " << V << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR NORM: " << V.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR: " << Vexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR NORM: " << Vexp.norm() << std::endl;

    // return success or failure based on the error
    return (V - Vexp).norm() > 1e-8;
}
