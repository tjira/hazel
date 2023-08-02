#include "../include/integral.h"

int test_integral_overlap_methane_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.21905884826794, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.18426110292136, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.07910270090866, 0.00000000000000, 0.00000000000000, 0.07910270090866, 0.00000000000000, 0.07910270090866, 0.03375905314694, 0.09298456485544, 0.03375901894552, 0.09298453890021, 0.03375898007309, 0.09298450940017, 0.03375896356095, 0.09298449686920, 0.21905884826794, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.71023506784998, -0.00000000000000, 0.00000000000000, 0.71023506784998, 0.00000000000000, 0.71023506784998, 0.27811166388164, 0.47930277453754, 0.27811151654629, 0.47930266970336, 0.27811134908876, 0.47930255055152, 0.27811127795650, 0.47930249993843, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.18184827232757, 0.12935163865812, 0.07054679336094, 0.05018110245515, 0.11503719063965, 0.08182789413294, -0.36743452700103, -0.26136242475621, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.05537966624479, 0.03939245880873, 0.28336565040227, 0.20156267999816, -0.32028581006577, -0.22782469923524, -0.01846164271043, -0.01313208027336, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.32183876061644, 0.22892915361024, -0.23332989002539, -0.16597141499130, -0.15459333467970, -0.10996484661603, 0.06608053811153, 0.04700421000440, 0.18426110292136, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.62993628835041, 0.37983546655708, 0.70768809236746, 0.37983537384829, 0.70768798779186, 0.37983526847786, 0.70768786893390, 0.37983522371881, 0.70768781844565, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.25079964621995, 0.28290222779780, 0.09729604249047, 0.10975003459845, 0.15865594263200, 0.17896408394025, -0.50675508603878, -0.57162034086587, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.07637796347580, 0.08615441188865, 0.39080949038779, 0.44083350148248, -0.44172886025022, -0.49827065733061, -0.02546176434888, -0.02872090052401, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.44387030060462, 0.50068610058644, -0.32180165550628, -0.36299259375438, -0.21321062433681, -0.24050182702120, 0.09113636938126, 0.10280193325383, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.33333333333333, 0.20079181035452, 0.37284077274898, 0.13211878197322, 0.35671360599498, 0.15230288307509, 0.36145354956910, 0.44999071627648, 0.43136230438949, 0.00000000000000, -0.00000000000000, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.02461860619857, 0.00578140744455, 0.04886845864517, 0.01147622274721, -0.09007004185606, -0.02115197542647, 0.01658268330390, 0.00389426478708, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.14307095445517, 0.03359863164124, -0.04023942939168, -0.00944978965415, -0.04347438346522, -0.01020948887985, -0.05935509928568, -0.01393890655969, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.33333333333333, 0.12744984986396, 0.35561722430931, 0.31624265371136, 0.39995308241243, 0.37072482201417, 0.41274758019056, 0.12078554232870, 0.35405201068519, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.04357050856539, 0.01023205215402, -0.16162991311384, -0.03795701638503, 0.12104110025501, 0.02842519360964, -0.00298227998604, -0.00070035637309, 0.07910270090866, 0.71023506784998, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.62993628835041, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 0.00000000000000, 0.33333333333333, 0.00000000000000, 1.00000000000000, 0.37316243154241, 0.41332010714022, 0.25304231927957, 0.38511117546988, 0.17837566708038, 0.36757646097466, 0.13062695096179, 0.35636315963453, 0.03375905314694, 0.27811166388164, 0.18184827232757, 0.05537966624479, 0.32183876061644, 0.37983546655708, 0.25079964621995, 0.07637796347580, 0.44387030060462, 0.20079181035452, 0.02461860619857, 0.14307095445517, 0.12744984986396, 0.04357050856539, 0.37316243154241, 1.00000000000000, 0.65829196968307, 0.01933295536092, 0.15201753068638, 0.01933329581569, 0.15201859852808, 0.01933316301108, 0.15201818198572, 0.09298456485544, 0.47930277453754, 0.12935163865812, 0.03939245880873, 0.22892915361024, 0.70768809236746, 0.28290222779780, 0.08615441188865, 0.50068610058644, 0.37284077274898, 0.00578140744455, 0.03359863164124, 0.35561722430931, 0.01023205215402, 0.41332010714022, 0.65829196968307, 1.00000000000000, 0.15201753068638, 0.40628653370048, 0.15201859852808, 0.40628829129637, 0.15201818198572, 0.40628760569613, 0.03375901894552, 0.27811151654629, 0.07054679336094, 0.28336565040227, -0.23332989002539, 0.37983537384829, 0.09729604249047, 0.39080949038779, -0.32180165550628, 0.13211878197322, 0.04886845864517, -0.04023942939168, 0.31624265371136, -0.16162991311384, 0.25304231927957, 0.01933295536092, 0.15201753068638, 1.00000000000000, 0.65829196968307, 0.01933361155280, 0.15201958883211, 0.01933358505535, 0.15201950572370, 0.09298453890021, 0.47930266970336, 0.05018110245515, 0.20156267999816, -0.16597141499130, 0.70768798779186, 0.10975003459845, 0.44083350148248, -0.36299259375438, 0.35671360599498, 0.01147622274721, -0.00944978965415, 0.39995308241243, -0.03795701638503, 0.38511117546988, 0.15201753068638, 0.40628653370048, 0.65829196968307, 1.00000000000000, 0.15201958883211, 0.40628992126619, 0.15201950572370, 0.40628978447583, 0.03375898007309, 0.27811134908876, 0.11503719063965, -0.32028581006577, -0.15459333467970, 0.37983526847786, 0.15865594263200, -0.44172886025022, -0.21321062433681, 0.15230288307509, -0.09007004185606, -0.04347438346522, 0.37072482201417, 0.12104110025501, 0.17837566708038, 0.01933329581569, 0.15201859852808, 0.01933361155280, 0.15201958883211, 1.00000000000000, 0.65829196968307, 0.01933357642783, 0.15201947866375, 0.09298450940017, 0.47930255055152, 0.08182789413294, -0.22782469923524, -0.10996484661603, 0.70768786893390, 0.17896408394025, -0.49827065733061, -0.24050182702120, 0.36145354956910, -0.02115197542647, -0.01020948887985, 0.41274758019056, 0.02842519360964, 0.36757646097466, 0.15201859852808, 0.40628829129637, 0.15201958883211, 0.40628992126619, 0.65829196968307, 1.00000000000000, 0.15201947866375, 0.40628973993714, 0.03375896356095, 0.27811127795650, -0.36743452700103, -0.01846164271043, 0.06608053811153, 0.37983522371881, -0.50675508603878, -0.02546176434888, 0.09113636938126, 0.44999071627648, 0.01658268330390, -0.05935509928568, 0.12078554232870, -0.00298227998604, 0.13062695096179, 0.01933316301108, 0.15201818198572, 0.01933358505535, 0.15201950572370, 0.01933357642783, 0.15201947866375, 1.00000000000000, 0.65829196968307, 0.09298449686920, 0.47930249993843, -0.26136242475621, -0.01313208027336, 0.04700421000440, 0.70768781844565, -0.57162034086587, -0.02872090052401, 0.10280193325383, 0.43136230438949, 0.00389426478708, -0.01393890655969, 0.35405201068519, -0.00070035637309, 0.35636315963453, 0.15201818198572, 0.40628760569613, 0.15201950572370, 0.40628978447583, 0.15201947866375, 0.40628973993714, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}