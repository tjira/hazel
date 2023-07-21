#include "../include/integral.h"

int test_int_overlap_methane_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/methane.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix S = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Sexp(system.shells.nbf(), system.shells.nbf()); Sexp << 1.00000000000000, 0.21905884826794, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.18426110292136, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.03375905314694, 0.09298456485544, 0.03375901894552, 0.09298453890021, 0.03375898007309, 0.09298450940017, 0.03375896356095, 0.09298449686920, 0.21905884826794, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.27811166388164, 0.47930277453754, 0.27811151654629, 0.47930266970336, 0.27811134908876, 0.47930255055152, 0.27811127795650, 0.47930249993843, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.18184827232757, 0.12935163865812, 0.07054679336094, 0.05018110245515, 0.11503719063965, 0.08182789413294, -0.36743452700103, -0.26136242475621, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.05537966624479, 0.03939245880873, 0.28336565040227, 0.20156267999816, -0.32028581006577, -0.22782469923524, -0.01846164271043, -0.01313208027336, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.32183876061644, 0.22892915361024, -0.23332989002539, -0.16597141499130, -0.15459333467970, -0.10996484661603, 0.06608053811153, 0.04700421000440, 0.18426110292136, 0.81227322234085, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.37983546655708, 0.70768809236746, 0.37983537384829, 0.70768798779186, 0.37983526847786, 0.70768786893390, 0.37983522371881, 0.70768781844565, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, -0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.25079964621995, 0.28290222779780, 0.09729604249047, 0.10975003459845, 0.15865594263200, 0.17896408394025, -0.50675508603878, -0.57162034086587, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.07637796347580, 0.08615441188865, 0.39080949038779, 0.44083350148248, -0.44172886025022, -0.49827065733061, -0.02546176434888, -0.02872090052401, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.56975428958387, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.44387030060462, 0.50068610058644, -0.32180165550628, -0.36299259375438, -0.21321062433681, -0.24050182702120, 0.09113636938126, 0.10280193325383, 0.03375905314694, 0.27811166388164, 0.18184827232757, 0.05537966624479, 0.32183876061644, 0.37983546655708, 0.25079964621995, 0.07637796347580, 0.44387030060462, 1.00000000000000, 0.65829196968307, 0.01933295536092, 0.15201753068638, 0.01933329581569, 0.15201859852808, 0.01933316301108, 0.15201818198572, 0.09298456485544, 0.47930277453754, 0.12935163865812, 0.03939245880873, 0.22892915361024, 0.70768809236746, 0.28290222779780, 0.08615441188865, 0.50068610058644, 0.65829196968307, 1.00000000000000, 0.15201753068638, 0.40628653370048, 0.15201859852808, 0.40628829129637, 0.15201818198572, 0.40628760569613, 0.03375901894552, 0.27811151654629, 0.07054679336094, 0.28336565040227, -0.23332989002539, 0.37983537384829, 0.09729604249047, 0.39080949038779, -0.32180165550628, 0.01933295536092, 0.15201753068638, 1.00000000000000, 0.65829196968307, 0.01933361155280, 0.15201958883211, 0.01933358505535, 0.15201950572370, 0.09298453890021, 0.47930266970336, 0.05018110245515, 0.20156267999816, -0.16597141499130, 0.70768798779186, 0.10975003459845, 0.44083350148248, -0.36299259375438, 0.15201753068638, 0.40628653370048, 0.65829196968307, 1.00000000000000, 0.15201958883211, 0.40628992126619, 0.15201950572370, 0.40628978447583, 0.03375898007309, 0.27811134908876, 0.11503719063965, -0.32028581006577, -0.15459333467970, 0.37983526847786, 0.15865594263200, -0.44172886025022, -0.21321062433681, 0.01933329581569, 0.15201859852808, 0.01933361155280, 0.15201958883211, 1.00000000000000, 0.65829196968307, 0.01933357642783, 0.15201947866375, 0.09298450940017, 0.47930255055152, 0.08182789413294, -0.22782469923524, -0.10996484661603, 0.70768786893390, 0.17896408394025, -0.49827065733061, -0.24050182702120, 0.15201859852808, 0.40628829129637, 0.15201958883211, 0.40628992126619, 0.65829196968307, 1.00000000000000, 0.15201947866375, 0.40628973993714, 0.03375896356095, 0.27811127795650, -0.36743452700103, -0.01846164271043, 0.06608053811153, 0.37983522371881, -0.50675508603878, -0.02546176434888, 0.09113636938126, 0.01933316301108, 0.15201818198572, 0.01933358505535, 0.15201950572370, 0.01933357642783, 0.15201947866375, 1.00000000000000, 0.65829196968307, 0.09298449686920, 0.47930249993843, -0.26136242475621, -0.01313208027336, 0.04700421000440, 0.70768781844565, -0.57162034086587, -0.02872090052401, 0.10280193325383, 0.15201818198572, 0.40628760569613, 0.15201950572370, 0.40628978447583, 0.15201947866375, 0.40628973993714, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP: " << S << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED OVERLAP NORM: " << S.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP: " << Sexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED OVERLAP NORM: " << Sexp.norm() << std::endl;

    // return success or failure based on the error
    return (S - Sexp).norm() > 1e-8;
}
