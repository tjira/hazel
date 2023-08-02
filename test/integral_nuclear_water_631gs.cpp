#include "../include/integral.h"

int test_integral_nuclear_water_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/water.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << -62.55089151305392, -5.17381526318127, -0.02362660844608, -0.00818531992164, -0.00944968693583, -5.49733846914846, -0.00325570439474, -0.00112792244820, -0.00130214996181, -0.63382736482973, -0.00003309024421, -0.00033398422854, -0.63421644700227, 0.00029412829430, -0.63339214029141, -0.96522242250009, -2.16096963942792, -0.96521570192472, -2.16096574949510, -5.17381526318127, -11.60063670115483, -0.16084613529901, -0.05572434136776, -0.06433194136088, -7.84446630757583, -0.10369738082441, -0.03592546530106, -0.04147474637218, -4.89492899880369, -0.00226830856100, -0.02289421614217, -4.92160013274041, 0.02016214165872, -4.86509487495136, -2.15098223209886, -3.61205963483215, -2.15097105965871, -3.61205373217492, -0.02362660844608, -0.16084613529901, -12.70649947735247, -0.00258042096741, -0.02604446736896, -0.15277789148369, -4.30932599931437, -0.00240249573412, -0.02424834991993, -0.19859634253308, -0.01419341821872, -0.03847083356817, -0.09193509135063, 0.01289322105441, -0.06834004431421, -1.64157342481444, -0.92957757283930, -1.01828317246575, -0.59883463624806, -0.00818531992164, -0.05572434136776, -0.00258042096741, -12.73684055510595, 0.02293645750952, -0.05292916882668, -0.00240249573412, -4.33757469923098, 0.02135468351697, -0.01419341821872, -0.09193509135063, 0.01289322105441, -0.09869891595233, -0.02881944965650, -0.01143695410548, 1.15821118374787, 0.59433474920633, -2.07970369361967, -1.12384839716832, -0.00944968693583, -0.06433194136088, -0.02604446736896, 0.02293645750952, -12.67256017287372, -0.06110496177916, -0.02424834991993, 0.02135468351697, -4.27772725195728, -0.03847083356817, 0.01289322105441, -0.06834004431421, -0.02881944965650, -0.01143695410548, -0.07624376139636, -1.15508547067890, -0.63633119360333, 0.09125008572557, 0.02502902220534, -5.49733846914846, -7.84446630757583, -0.15277789148369, -0.05292916882668, -0.06110496177916, -7.64761238060573, -0.22302394807257, -0.07726582283723, -0.08920046561539, -5.08345078632496, -0.00400171985532, -0.04038896678111, -5.13050294033014, 0.03556917143373, -5.03081870285366, -2.87283274249483, -4.41771876105441, -2.87282376339514, -4.41771296491345, -0.00325570439474, -0.10369738082441, -4.30932599931437, -0.00240249573412, -0.02424834991993, -0.22302394807257, -5.35543935935863, -0.00583186671613, -0.05885878572933, -0.25560232875281, -0.01360017362047, -0.05591826316201, -0.13206742133351, 0.02347892223782, -0.08910006876437, -2.45581726448258, -1.78191546750052, -1.52145881595374, -1.14628226301400, -0.00112792244820, -0.03592546530106, -0.00240249573412, -4.33757469923098, 0.02135468351697, -0.07726582283723, -0.00583186671613, -5.42400867244711, 0.05183492834965, -0.01360017362047, -0.13206742133351, 0.02347892223782, -0.14299424687878, -0.03834281396456, -0.00858052247551, 1.73799704645614, 1.14380076868244, -3.11590962979468, -2.15826702233689, -0.00130214996181, -0.04147474637218, -0.02424834991993, 0.02135468351697, -4.27772725195728, -0.08920046561539, -0.05885878572933, 0.05183492834965, -5.27873860760021, -0.05591826316201, 0.02347892223782, -0.08910006876437, -0.03834281396456, -0.00858052247551, -0.09642742530911, -1.72955334423892, -1.22109233137843, 0.13880543289841, 0.04993375343181, -0.63382736482973, -4.89492899880369, -0.19859634253308, -0.01419341821872, -0.03847083356817, -5.08345078632496, -0.25560232875281, -0.01360017362047, -0.05591826316201, -7.13495631924680, 0.00876837333204, -0.04551924825517, -2.40918604096685, 0.02245311927718, -2.36996642638257, -2.44467861966937, -2.93216593076883, -1.62220181188359, -2.77069149611053, -0.00003309024421, -0.00226830856100, -0.01419341821872, -0.09193509135063, 0.01289322105441, -0.00400171985532, -0.01360017362047, -0.13206742133351, 0.02347892223782, 0.00876837333204, -2.40918604096685, 0.02245311927718, -0.02736913122681, -0.01716948814386, 0.00878571015258, 0.94507860801817, 0.17737517745734, -1.04878379049335, -0.21406768869385, -0.00033398422854, -0.02289421614217, -0.03847083356817, 0.01289322105441, -0.06834004431421, -0.04038896678111, -0.05591826316201, 0.02347892223782, -0.08910006876437, -0.04551924825517, 0.02245311927718, -2.36996642638257, -0.01716948814386, 0.00878571015258, -0.03637417947363, -0.94330633219498, -0.20726230369988, 0.04348557537435, -0.01353172796402, -0.63421644700227, -4.92160013274041, -0.09193509135063, -0.09869891595233, -0.02881944965650, -5.13050294033014, -0.13206742133351, -0.14299424687878, -0.03834281396456, -2.40918604096685, -0.02736913122681, -0.01716948814386, -7.24632716386093, 0.03427517316559, -2.37400134685177, -1.78553845193563, -2.80224112570359, -3.26578352934576, -3.09284626017656, 0.00029412829430, 0.02016214165872, 0.01289322105441, -0.02881944965650, -0.01143695410548, 0.03556917143373, 0.02347892223782, -0.03834281396456, -0.00858052247551, 0.02245311927718, -0.01716948814386, 0.00878571015258, 0.03427517316559, -2.37400134685177, 0.03051299582406, 0.66936201249611, 0.14128212033882, 0.09980378032484, 0.02946443488870, -0.63339214029141, -4.86509487495136, -0.06834004431421, -0.01143695410548, -0.07624376139636, -5.03081870285366, -0.08910006876437, -0.00858052247551, -0.09642742530911, -2.36996642638257, 0.00878571015258, -0.03637417947363, -2.37400134685177, 0.03051299582406, -7.04104916511632, -1.76848264877286, -2.76444029499312, -1.11069334413613, -2.63529760929621, -0.96522242250009, -2.15098223209886, -1.64157342481444, 1.15821118374787, -1.15508547067890, -2.87283274249483, -2.45581726448258, 1.73799704645614, -1.72955334423892, -2.44467861966937, 0.94507860801817, -0.94330633219498, -1.78553845193563, 0.66936201249611, -1.76848264877286, -6.27313096116240, -3.74384427038999, -0.39369830641728, -1.45518634103685, -2.16096963942792, -3.61205963483215, -0.92957757283930, 0.59433474920633, -0.63633119360333, -4.41771876105441, -1.78191546750052, 1.14380076868244, -1.22109233137843, -2.93216593076883, 0.17737517745734, -0.20726230369988, -2.80224112570359, 0.14128212033882, -2.76444029499312, -3.74384427038999, -4.69136306286849, -1.45518502976856, -2.82583716583917, -0.96521570192472, -2.15097105965871, -1.01828317246575, -2.07970369361967, 0.09125008572557, -2.87282376339514, -1.52145881595374, -3.11590962979468, 0.13880543289841, -1.62220181188359, -1.04878379049335, 0.04348557537435, -3.26578352934576, 0.09980378032484, -1.11069334413613, -0.39369830641728, -1.45518502976856, -6.27312419839383, -3.74384026304667, -2.16096574949510, -3.61205373217492, -0.59883463624806, -1.12384839716832, 0.02502902220534, -4.41771296491345, -1.14628226301400, -2.15826702233689, 0.04993375343181, -2.77069149611053, -0.21406768869385, -0.01353172796402, -3.09284626017656, 0.02946443488870, -2.63529760929621, -1.45518634103685, -2.82583716583917, -3.74384026304667, -4.69135975903299;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}