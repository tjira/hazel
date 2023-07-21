#include "../include/integral.h"

int test_int_overlap_formaldehyde_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix S = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Sexp(system.shells.nbf(), system.shells.nbf()); Sexp << 1.00000000000000, 0.23368985719701, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.16727976258450, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.03353153616877, 0.00000000000000, 0.00000000000000, 0.03353153616877, 0.00000000000000, 0.03353153616877, 0.00000250692049, 0.02069370383582, -0.04291572176216, -0.00108572491034, -0.00494133758277, 0.04992224315450, -0.09178187650130, -0.00232198983365, -0.01056781098281, 0.05897531479128, 0.00147261046013, 0.00670212624038, 0.00080432611406, 0.00016955710198, 0.00153875669385, 0.00003811784665, 0.01115302903619, 0.00003811805197, 0.01115304516041, 0.23368985719701, 1.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.76364080963424, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.54706575927831, 0.00000000000000, 0.00000000000000, 0.54706575927831, -0.00000000000000, 0.54706575927831, 0.00790400131029, 0.16995248977228, -0.26857693595286, -0.00679472831247, -0.03092408220103, 0.28684667448930, -0.45716095140195, -0.01156571560706, -0.05263773968565, 0.34441626279181, 0.00737484049599, 0.03356428148850, 0.05309566267203, 0.00084914280859, 0.05677369466218, 0.00264321356744, 0.07757228178486, 0.00264322302022, 0.07757237971530, 0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.01883490737296, 0.20029043145586, -0.27049618124563, -0.00908250348482, -0.04133617584686, 0.14539926353993, -0.08168666334969, -0.00598850726125, -0.02725481907328, 0.36224555820974, 0.00985722244979, 0.04486206707513, 0.07263343804951, 0.00142539607234, 0.07880749063484, 0.00440368681109, 0.05187797252603, 0.00492129818897, 0.05797559760252, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00047650435037, 0.00506714793101, -0.00908250348482, 0.08828044141164, -0.00104576397582, 0.00367845619015, -0.00598850726125, 0.15487106696564, -0.00068951970932, 0.01420925399365, -0.09939040800694, 0.00142539607234, -0.00320724960375, -0.01144386109066, 0.00199375082582, -0.00032355279884, -0.00381163873114, 0.00055946754814, 0.00659083522214, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50152068502883, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00216866061780, 0.02306154005535, -0.04133617584686, -0.00104576397582, 0.08375075180195, 0.01674134363670, -0.02725481907328, -0.00068951970932, 0.15188443677450, 0.06466897840619, 0.00142539607234, -0.09321635542161, 0.00836305024040, -0.00235828071804, -0.01388591064538, 0.00288152211260, 0.03394599375548, -0.00180784906938, -0.02129745570130, 0.16727976258450, 0.76364080963424, 0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.00000000000000, 0.69901474023950, 0.00000000000000, 0.69901474023950, 0.06598002437718, 0.37105640729895, -0.33271618783189, -0.00841738734362, -0.03830910761427, 0.55432172289969, -0.63994782355550, -0.01619004096457, -0.07368391120329, 0.41126894069658, 0.00432190390863, 0.01966979481583, 0.24054529711307, 0.00049762616906, 0.24270074723338, 0.04547362136610, 0.21373753718753, 0.04547370168643, 0.21373773442623, -0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.14780640398267, 0.57181730489478, -0.26727953671892, -0.01385604370389, -0.06306145216930, 0.50586447092026, -0.04466007052387, -0.01477475294289, -0.06724267012522, 0.44173373083636, 0.00402559730873, 0.01832124793792, 0.42678424044988, 0.00088290650093, 0.43060851867714, 0.10866290088618, 0.26683829509837, 0.12143496356453, 0.29820178523759, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00373935443988, 0.01446640686902, -0.01385604370389, 0.28006124567440, -0.00159539177465, 0.01279786602167, -0.01477475294289, 0.53897152073324, -0.00170117241410, 0.01846038684042, -0.14378331799536, 0.00088290650093, 0.00351225064035, -0.01655528286169, 0.01089396557104, -0.00798380702807, -0.01960545354846, 0.01380508124379, 0.03390044968441, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.50152068502883, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.01701850298643, 0.06583932934451, -0.06306145216930, -0.00159539177465, 0.27315086272747, 0.05824548718534, -0.06724267012522, -0.00170117241410, 0.53160295269438, 0.08401668085365, 0.00088290650093, -0.13995903976809, 0.04914015005404, -0.00354082396064, 0.01642524032504, 0.07110282024002, 0.17460379922496, -0.04460938505258, -0.10954504263479, 0.03353153616877, 0.54706575927831, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.69901474023950, 0.00000000000000, -0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.33333333333333, 0.09512651348186, 0.40892000124827, -0.44563295250026, -0.01860307834326, -0.08466609663016, 0.39399043541800, -0.41329324047328, -0.01544118143198, -0.07027571109904, 0.48500800352955, 0.01346222033787, 0.06126908820958, 0.20857444070816, 0.00252759111166, 0.21952261205538, 0.01865976770461, 0.13283389347137, 0.02263472815746, 0.14468136553304, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00234474668007, 0.00692209847914, -0.01015025852770, 0.14449830367258, -0.00159063873146, 0.00233790256671, -0.00112910748015, 0.09843548540738, -0.00041700952761, 0.01346222033787, -0.12722723864429, 0.00155004675610, -0.00423367330302, -0.01926824079112, 0.00447975725504, -0.00117341095324, -0.00349739447713, 0.00226746926390, 0.00675826159993, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.01067138165656, 0.03150376770469, -0.04619572919439, -0.00159063873146, 0.13760850836385, 0.01064023281376, -0.00513877978992, -0.00041700952761, 0.09662922341326, 0.06126908820958, 0.00155004675610, -0.12051326603864, 0.01937140095776, -0.00401025599961, -0.01825142676443, 0.01045025609742, 0.03114737241788, -0.00732704195665, -0.02183847299944, 0.03353153616877, 0.54706575927831, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.69901474023950, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.33333333333333, 0.00250444113510, 0.13548359484381, -0.18961843239452, 0.00253185258753, -0.02183276076702, 0.30163871891303, -0.46728146343033, -0.00683647934567, -0.05380302048228, 0.20857444070816, -0.00423367330302, 0.01937140095776, 0.12046237547328, -0.00048746725317, 0.04256345305548, 0.00277537570778, 0.08548987926303, 0.00294694277105, 0.08600131809728, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00026997530081, 0.00079701386723, -0.00159063873146, 0.01663760666875, 0.00348135786306, 0.00026918726619, -0.00041700952761, 0.01133391082684, 0.00244462287057, 0.00252759111166, -0.01926824079112, -0.00401025599961, -0.00048746725317, 0.03789967659734, -0.00046174287356, -0.00076781336956, -0.00228849597045, -0.00083295952433, -0.00248266137814, 0.03353153616877, 0.54706575927831, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.69901474023950, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 1.00000000000000, 0.00367382958530, 0.13893583217139, -0.19650822770325, -0.00497146195250, 0.01072964124976, 0.30280469401919, -0.46908772542445, -0.01186745108125, -0.03132207288506, 0.21952261205538, 0.00447975725504, -0.01825142676443, 0.04256345305548, -0.00046174287356, 0.11635090543998, 0.00952721454291, 0.10561398308814, 0.00538077426290, 0.09325542741375, 0.00000250692049, 0.00790400131029, 0.01883490737296, 0.00047650435037, 0.00216866061780, 0.06598002437718, 0.14780640398267, 0.00373935443988, 0.01701850298643, 0.09512651348186, 0.00234474668007, 0.01067138165656, 0.00250444113510, 0.00026997530081, 0.00367382958530, 1.00000000000000, 0.21905884826794, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.18426110292136, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.07910270090866, 0.00000000000000, 0.00000000000000, 0.07910270090866, 0.00000000000000, 0.07910270090866, 0.03113396272593, 0.09093209448755, 0.03113396152718, 0.09093209352146, 0.02069370383582, 0.16995248977228, 0.20029043145586, 0.00506714793101, 0.02306154005535, 0.37105640729895, 0.57181730489478, 0.01446640686902, 0.06583932934451, 0.40892000124827, 0.00692209847914, 0.03150376770469, 0.13548359484381, 0.00079701386723, 0.13893583217139, 0.21905884826794, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.71023506784998, -0.00000000000000, 0.00000000000000, 0.71023506784998, 0.00000000000000, 0.71023506784998, 0.26657345989387, 0.47099242546923, 0.26657345451603, 0.47099242154791, -0.04291572176216, -0.26857693595286, -0.27049618124563, -0.00908250348482, -0.04133617584686, -0.33271618783189, -0.26727953671892, -0.01385604370389, -0.06306145216930, -0.44563295250026, -0.01015025852770, -0.04619572919439, -0.18961843239452, -0.00159063873146, -0.19650822770325, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.16164573131715, 0.11851147137673, 0.22734329118236, 0.16667800697713, -0.00108572491034, -0.00679472831247, -0.00908250348482, 0.08828044141164, -0.00104576397582, -0.00841738734362, -0.01385604370389, 0.28006124567440, -0.00159539177465, -0.01860307834326, 0.14449830367258, -0.00159063873146, 0.00253185258753, 0.01663760666875, -0.00497146195250, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.05111972102527, -0.03747870918525, 0.06096079962707, 0.04469375160678, -0.00494133758277, -0.03092408220103, -0.04133617584686, -0.00104576397582, 0.08375075180195, -0.03830910761427, -0.06306145216930, -0.00159539177465, 0.27315086272747, -0.08466609663016, -0.00159063873146, 0.13760850836385, -0.02183276076702, 0.00348135786306, 0.01072964124976, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.32000197225035, 0.23461123449302, -0.27521379954899, -0.20177453824504, 0.04992224315450, 0.28684667448930, 0.14539926353993, 0.00367845619015, 0.01674134363670, 0.55432172289969, 0.50586447092026, 0.01279786602167, 0.05824548718534, 0.39399043541800, 0.00233790256671, 0.01064023281376, 0.30163871891303, 0.00026918726619, 0.30280469401919, 0.18426110292136, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.62993628835041, 0.37249426861346, 0.69937508117444, 0.37249426515321, 0.69937507724084, -0.09178187650130, -0.45716095140195, -0.08168666334969, -0.00598850726125, -0.02725481907328, -0.63994782355550, -0.04466007052387, -0.01477475294289, -0.06724267012522, -0.41329324047328, -0.00112910748015, -0.00513877978992, -0.46728146343033, -0.00041700952761, -0.46908772542445, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.22947299077555, 0.26086373498577, 0.32273753999303, 0.36688640405716, -0.00232198983365, -0.01156571560706, -0.00598850726125, 0.15487106696564, -0.00068951970932, -0.01619004096457, -0.01477475294289, 0.53897152073324, -0.00170117241410, -0.01544118143198, 0.09843548540738, -0.00041700952761, -0.00683647934567, 0.01133391082684, -0.01186745108125, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.07256978069075, -0.08249695955113, 0.08654022032199, 0.09837848500963, -0.01056781098281, -0.05263773968565, -0.02725481907328, -0.00068951970932, 0.15188443677450, -0.07368391120329, -0.06724267012522, -0.00170117241410, 0.53160295269438, -0.07027571109904, -0.00041700952761, 0.09662922341326, -0.05380302048228, 0.00244462287057, -0.03132207288506, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.45427620654139, 0.51641889336539, -0.39069472504173, -0.44413978850356, 0.05897531479128, 0.34441626279181, 0.36224555820974, 0.01420925399365, 0.06466897840619, 0.41126894069658, 0.44173373083636, 0.01846038684042, 0.08401668085365, 0.48500800352955, 0.01346222033787, 0.06126908820958, 0.20857444070816, 0.00252759111166, 0.21952261205538, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.33333333333333, 0.18027758347407, 0.36332825851084, 0.24522357140857, 0.37918505125145, 0.00147261046013, 0.00737484049599, 0.00985722244979, -0.09939040800694, 0.00142539607234, 0.00432190390863, 0.00402559730873, -0.14378331799536, 0.00088290650093, 0.01346222033787, -0.12722723864429, 0.00155004675610, -0.00423367330302, -0.01926824079112, 0.00447975725504, 0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.02099995967951, -0.00512721506972, 0.03522074607333, 0.00859927095688, 0.00670212624038, 0.03356428148850, 0.04486206707513, 0.00142539607234, -0.09321635542161, 0.01966979481583, 0.01832124793792, 0.00088290650093, -0.13995903976809, 0.06126908820958, 0.00155004675610, -0.12051326603864, 0.01937140095776, -0.00401025599961, -0.01825142676443, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.13145667425098, 0.03209561596887, -0.15900768049453, -0.03882229314366, 0.00080432611406, 0.05309566267203, 0.07263343804951, -0.00320724960375, 0.00836305024040, 0.24054529711307, 0.42678424044988, 0.00351225064035, 0.04914015005404, 0.20857444070816, -0.00423367330302, 0.01937140095776, 0.12046237547328, -0.00048746725317, 0.04256345305548, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.33333333333333, 0.12051472663978, 0.34873694382039, 0.12331782395927, 0.34942132758317, 0.00016955710198, 0.00084914280859, 0.00142539607234, -0.01144386109066, -0.00235828071804, 0.00049762616906, 0.00088290650093, -0.01655528286169, -0.00354082396064, 0.00252759111166, -0.01926824079112, -0.00401025599961, -0.00048746725317, 0.03789967659734, -0.00046174287356, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, -0.04157257021180, -0.01015009132065, -0.04263699755281, -0.01040997524530, 0.00153875669385, 0.05677369466218, 0.07880749063484, 0.00199375082582, -0.01388591064538, 0.24270074723338, 0.43060851867714, 0.01089396557104, 0.01642524032504, 0.21952261205538, 0.00447975725504, -0.01825142676443, 0.04256345305548, -0.00046174287356, 0.11635090543998, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 1.00000000000000, 0.37411179110540, 0.41065357177632, 0.30636269344212, 0.39411238627546, 0.00003811784665, 0.00264321356744, 0.00440368681109, -0.00032355279884, 0.00288152211260, 0.04547362136610, 0.10866290088618, -0.00798380702807, 0.07110282024002, 0.01865976770461, -0.00117341095324, 0.01045025609742, 0.00277537570778, -0.00076781336956, 0.00952721454291, 0.03113396272593, 0.26657345989387, 0.16164573131715, -0.05111972102527, 0.32000197225035, 0.37249426861346, 0.22947299077555, -0.07256978069075, 0.45427620654139, 0.18027758347407, -0.02099995967951, 0.13145667425098, 0.12051472663978, -0.04157257021180, 0.37411179110540, 1.00000000000000, 0.65829196968307, 0.01352826408506, 0.13179550057842, 0.01115302903619, 0.07757228178486, 0.05187797252603, -0.00381163873114, 0.03394599375548, 0.21373753718753, 0.26683829509837, -0.01960545354846, 0.17460379922496, 0.13283389347137, -0.00349739447713, 0.03114737241788, 0.08548987926303, -0.00228849597045, 0.10561398308814, 0.09093209448755, 0.47099242546923, 0.11851147137673, -0.03747870918525, 0.23461123449302, 0.69937508117444, 0.26086373498577, -0.08249695955113, 0.51641889336539, 0.36332825851084, -0.00512721506972, 0.03209561596887, 0.34873694382039, -0.01015009132065, 0.41065357177632, 0.65829196968307, 1.00000000000000, 0.13179550057842, 0.37208891461058, 0.00003811805197, 0.00264322302022, 0.00492129818897, 0.00055946754814, -0.00180784906938, 0.04547370168643, 0.12143496356453, 0.01380508124379, -0.04460938505258, 0.02263472815746, 0.00226746926390, -0.00732704195665, 0.00294694277105, -0.00083295952433, 0.00538077426290, 0.03113396152718, 0.26657345451603, 0.22734329118236, 0.06096079962707, -0.27521379954899, 0.37249426515321, 0.32273753999303, 0.08654022032199, -0.39069472504173, 0.24522357140857, 0.03522074607333, -0.15900768049453, 0.12331782395927, -0.04263699755281, 0.30636269344212, 0.01352826408506, 0.13179550057842, 1.00000000000000, 0.65829196968307, 0.01115304516041, 0.07757237971530, 0.05797559760252, 0.00659083522214, -0.02129745570130, 0.21373773442623, 0.29820178523759, 0.03390044968441, -0.10954504263479, 0.14468136553304, 0.00675826159993, -0.02183847299944, 0.08600131809728, -0.00248266137814, 0.09325542741375, 0.09093209352146, 0.47099242154791, 0.16667800697713, 0.04469375160678, -0.20177453824504, 0.69937507724084, 0.36688640405716, 0.09837848500963, -0.44413978850356, 0.37918505125145, 0.00859927095688, -0.03882229314366, 0.34942132758317, -0.01040997524530, 0.39411238627546, 0.13179550057842, 0.37208891461058, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP: " << S << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP NORM: " << S.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP: " << Sexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP NORM: " << Sexp.norm() << std::endl;

    // return success or failure based on the error
    return (S - Sexp).norm() > 1e-8;
}
