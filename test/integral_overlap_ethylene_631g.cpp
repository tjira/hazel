#include "../include/integral.h"

int test_integral_overlap_ethylene_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ethylene.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.21905884826794, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.18426110292136, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000874324188, 0.02335378654697, -0.04810845278834, -0.00066743946116, -0.01274306189421, 0.06840394817287, -0.12924290077857, -0.00179306976331, -0.03423411455893, 0.03391670304205, 0.09310399695649, 0.03391651020287, 0.09310385111834, 0.00005903151461, 0.01561843030530, 0.00005902641554, 0.01561804698325, 0.21905884826794, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81227322234085, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.02335378654697, 0.22365959146726, -0.29338985806257, -0.00407038591817, -0.07771368447192, 0.37675917134285, -0.55304502475544, -0.00767274879829, -0.14649165733413, 0.27879000314479, 0.47978509720995, 0.27878917436284, 0.47978450833049, 0.00935998498366, 0.11792630370705, 0.00935954932402, 0.11792402231485, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.04810845278834, 0.29338985806257, -0.31552861778788, -0.00617585120993, -0.11791219848009, 0.24986148137832, -0.12740001061654, -0.00530540214780, -0.10129318369290, -0.17646880632742, -0.12530631990004, -0.20900679041200, -0.14841108781330, 0.01768584845567, 0.10701128289311, 0.01679449959120, 0.10162081316276, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00066743946116, 0.00407038591817, -0.00617585120993, 0.12953565605273, -0.00163587165365, 0.00346648879417, -0.00530540214780, 0.25493513219650, -0.00140530538865, -0.31597842393792, -0.22436879523066, 0.31064490758392, 0.22058206131634, -0.00833925915723, -0.05045812888369, 0.00881764308350, 0.05335413868441, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.01274306189421, 0.07771368447192, -0.11791219848009, -0.00163587165365, 0.09838852182450, 0.06618380217281, -0.10129318369290, -0.00140530538865, 0.22817799909615, -0.09611791691236, -0.06825105635043, -0.00612715957044, -0.00435076016077, 0.00333276050445, 0.02016544346464, 0.00579658284247, 0.03507418954750, 0.18426110292136, 0.81227322234085, -0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.06840394817287, 0.37675917134285, -0.24986148137832, -0.00346648879417, -0.06618380217281, 0.59822353178206, -0.58614861334744, -0.00813201614221, -0.15526020119490, 0.38026202922266, 0.70816914798195, 0.38026150838943, 0.70816856074378, 0.08009231497702, 0.27515247541061, 0.08009060324170, 0.27514889845030, -0.00000000000000, -0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.12924290077857, 0.55304502475544, -0.12740001061654, -0.00530540214780, -0.10129318369290, 0.58614861334744, 0.02390610915299, -0.00796787443529, -0.15212633205407, -0.24297496914836, -0.27395359807989, -0.28777619604051, -0.32446703293565, 0.18644595731070, 0.38960870387733, 0.17705387261722, 0.36998555654130, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00179306976331, 0.00767274879829, -0.00530540214780, 0.25493513219650, -0.00140530538865, 0.00813201614221, -0.00796787443529, 0.59811298834299, -0.00211054630131, -0.43506186393881, -0.49053103466228, 0.42771916475843, 0.48225242472571, -0.08791329184628, -0.18370891052756, 0.09295887899560, 0.19425411075306, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.03423411455893, 0.14649165733413, -0.10129318369290, -0.00140530538865, 0.22817799909615, 0.15526020119491, -0.15212633205407, -0.00211054630131, 0.55792800879708, -0.13234207440070, -0.14921531870754, -0.00843633199782, -0.00951194591441, 0.03513428967219, 0.07341872818417, 0.06110973623434, 0.12769966246170, 0.00000874324188, 0.02335378654697, 0.04810845278834, 0.00066743946116, 0.01274306189421, 0.06840394817287, 0.12924290077857, 0.00179306976331, 0.03423411455893, 1.00000000000000, 0.21905884826794, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.18426110292136, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00005903336454, 0.01561856936807, 0.00005902612215, 0.01561802492663, 0.03391648716227, 0.09310383369343, 0.03391669526717, 0.09310399107661, 0.02335378654697, 0.22365959146726, 0.29338985806257, 0.00407038591817, 0.07771368447192, 0.37675917134285, 0.55304502475544, 0.00767274879829, 0.14649165733413, 0.21905884826794, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00936014303601, 0.11792713135468, 0.00935952425615, 0.11792389104153, 0.27878907533910, 0.47978443797048, 0.27878996973009, 0.47978507346762, -0.04810845278834, -0.29338985806257, -0.31552861778788, -0.00617585120993, -0.11791219848009, -0.24986148137832, -0.12740001061654, -0.00530540214780, -0.10129318369290, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.01679443532373, -0.10161662131102, -0.01768418796207, -0.10700434334394, 0.20903495559533, 0.14843112505171, 0.17651077310533, 0.12533613027607, -0.00066743946116, -0.00407038591817, -0.00617585120993, 0.12953565605273, -0.00163587165365, -0.00346648879417, -0.00530540214780, 0.25493513219650, -0.00140530538865, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, -0.00881778702559, -0.05335301293011, 0.00833928914948, 0.05045977577874, -0.31062840952521, -0.22057040253166, 0.31599393780091, 0.22437983052191, -0.01274306189421, -0.07771368447192, -0.11791219848009, -0.00163587165365, 0.09838852182450, -0.06618380217281, -0.10129318369290, -0.00140530538865, 0.22817799909615, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, -0.00580046945552, -0.03509639334253, -0.00333620326703, -0.02018686074903, 0.00599530884663, 0.00425713696833, 0.09598964138775, 0.06815997679048, 0.06840394817287, 0.37675917134285, 0.24986148137832, 0.00346648879417, 0.06618380217281, 0.59822353178206, 0.58614861334744, 0.00813201614221, 0.15526020119491, 0.18426110292136, 0.81227322234085, 0.00000000000000, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.08009293596438, 0.27515377306124, 0.08009050474705, 0.27514869262817, 0.38026144615967, 0.70816849057986, 0.38026200822380, 0.70816912430580, -0.12924290077857, -0.55304502475544, -0.12740001061654, -0.00530540214780, -0.10129318369290, -0.58614861334744, 0.02390610915299, -0.00796787443529, -0.15212633205407, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, -0.17704692244466, -0.36996678311249, -0.18643357782724, -0.38958630775133, 0.28781504595961, 0.32451085726777, 0.24303277197399, 0.27401877662571, -0.00179306976331, -0.00767274879829, -0.00530540214780, 0.25493513219650, -0.00140530538865, -0.00813201614221, -0.00796787443529, 0.59811298834299, -0.00211054630131, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.09295710308567, -0.19424816834537, 0.08791602509585, 0.18371625974482, -0.42769655299631, -0.48222696141736, 0.43508326024331, 0.49055516971577, -0.03423411455893, -0.14649165733413, -0.10129318369290, -0.00140530538865, 0.22817799909615, -0.15526020119490, -0.15212633205407, -0.00211054630131, 0.55792800879708, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.06114854390980, -0.12777929020139, -0.03517155058318, -0.07349724598585, 0.00825479205772, 0.00930726063429, 0.13216546625929, 0.14901619679676, 0.03391670304205, 0.27879000314479, -0.17646880632742, -0.31597842393792, -0.09611791691236, 0.38026202922266, -0.24297496914836, -0.43506186393881, -0.13234207440070, 0.00005903336454, 0.00936014303601, -0.01679443532373, -0.00881778702559, -0.00580046945552, 0.08009293596438, -0.17704692244466, -0.09295710308567, -0.06114854390980, 1.00000000000000, 0.65829196968307, 0.01482958008564, 0.13673308171036, 0.00066583606854, 0.03897894240830, 0.00001432990178, 0.00818403944996, 0.09310399695649, 0.47978509720995, -0.12530631990004, -0.22436879523066, -0.06825105635043, 0.70816914798195, -0.27395359807989, -0.49053103466228, -0.14921531870754, 0.01561856936807, 0.11792713135468, -0.10161662131102, -0.05335301293011, -0.03509639334253, 0.27515377306124, -0.36996678311249, -0.19424816834537, -0.12777929020139, 0.65829196968307, 1.00000000000000, 0.13673308171036, 0.38061594614700, 0.03897894240830, 0.17546788878124, 0.00818403944996, 0.06678427491659, 0.03391651020287, 0.27878917436284, -0.20900679041200, 0.31064490758392, -0.00612715957044, 0.38026150838943, -0.28777619604051, 0.42771916475843, -0.00843633199782, 0.00005902612215, 0.00935952425615, -0.01768418796207, 0.00833928914948, -0.00333620326703, 0.08009050474705, -0.18643357782724, 0.08791602509585, -0.03517155058318, 0.01482958008564, 0.13673308171036, 1.00000000000000, 0.65829196968307, 0.00001432950184, 0.00818394662591, 0.00066568220085, 0.03897528107064, 0.09310385111834, 0.47978450833049, -0.14841108781330, 0.22058206131634, -0.00435076016077, 0.70816856074378, -0.32446703293565, 0.48225242472571, -0.00951194591441, 0.01561802492663, 0.11792389104153, -0.10700434334394, 0.05045977577874, -0.02018686074903, 0.27514869262817, -0.38958630775133, 0.18371625974482, -0.07349724598585, 0.13673308171036, 0.38061594614700, 0.65829196968307, 1.00000000000000, 0.00818394662591, 0.06678380543495, 0.03897528107064, 0.17545770429861, 0.00005903151461, 0.00935998498366, 0.01768584845567, -0.00833925915723, 0.00333276050445, 0.08009231497702, 0.18644595731070, -0.08791329184628, 0.03513428967219, 0.03391648716227, 0.27878907533910, 0.20903495559533, -0.31062840952521, 0.00599530884663, 0.38026144615967, 0.28781504595961, -0.42769655299631, 0.00825479205772, 0.00066583606854, 0.03897894240830, 0.00001432950184, 0.00818394662591, 1.00000000000000, 0.65829196968307, 0.01482976357243, 0.13673375884475, 0.01561843030530, 0.11792630370705, 0.10701128289311, -0.05045812888369, 0.02016544346464, 0.27515247541061, 0.38960870387733, -0.18370891052756, 0.07341872818417, 0.09310383369343, 0.47978443797048, 0.14843112505171, -0.22057040253166, 0.00425713696833, 0.70816849057986, 0.32451085726777, -0.48222696141736, 0.00930726063429, 0.03897894240830, 0.17546788878124, 0.00818394662591, 0.06678380543495, 0.65829196968307, 1.00000000000000, 0.13673375884475, 0.38061710728974, 0.00005902641554, 0.00935954932402, 0.01679449959120, 0.00881764308350, 0.00579658284247, 0.08009060324170, 0.17705387261722, 0.09295887899560, 0.06110973623434, 0.03391669526717, 0.27878996973009, 0.17651077310533, 0.31599393780091, 0.09598964138775, 0.38026200822380, 0.24303277197399, 0.43508326024331, 0.13216546625929, 0.00001432990178, 0.00818403944996, 0.00066568220085, 0.03897528107064, 0.01482976357243, 0.13673375884475, 1.00000000000000, 0.65829196968307, 0.01561804698325, 0.11792402231485, 0.10162081316276, 0.05335413868441, 0.03507418954750, 0.27514889845030, 0.36998555654130, 0.19425411075306, 0.12769966246170, 0.09310399107661, 0.47978507346762, 0.12533613027607, 0.22437983052191, 0.06815997679048, 0.70816912430580, 0.27401877662571, 0.49055516971577, 0.14901619679676, 0.00818403944996, 0.06678427491659, 0.03897528107064, 0.17545770429861, 0.13673375884475, 0.38061710728974, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}