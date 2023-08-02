#include "../include/integral.h"

int test_integral_overlap_ammonia_631g(int, char**) {
    // initialize the system
    System system("../example/molecule/ammonia.xyz", "6-31G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Overlap(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 1.00000000000000, 0.22221581591139, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.17186916621874, -0.00000000000000, 0.00000000000000, -0.00000000000000, 0.03194056684349, 0.07782687765506, 0.03194059360971, 0.07782689516847, 0.03194058979707, 0.07782689267382, 0.22221581591139, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.77550678424825, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.24852804829932, 0.41673483339289, 0.24852817585408, 0.41673491190102, 0.24852815768492, 0.41673490071816, -0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.51929555839097, 0.00000000000000, 0.00000000000000, 0.10406281520907, 0.06432613555760, 0.26133836901302, 0.16154551380969, -0.14739774471992, -0.09111347019718, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.51929555839097, 0.00000000000000, 0.22752003215383, 0.14064086581738, -0.11147377517723, -0.06890717331451, 0.15535813444705, 0.09603416103637, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.51929555839097, 0.19164998940535, 0.11846790011719, -0.13638212535769, -0.08430419382570, -0.23121845462600, -0.14292698856849, 0.17186916621874, 0.77550678424825, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.39855150613730, 0.69577154167756, 0.39855159927097, 0.69577163288244, 0.39855158600482, 0.69577161989103, -0.00000000000000, 0.00000000000000, 0.51929555839097, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.00000000000000, 0.18212009644858, 0.17834840338088, 0.45736753043714, 0.44789543786962, -0.25796038297464, -0.25261802096777, 0.00000000000000, -0.00000000000000, 0.00000000000000, 0.51929555839097, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.00000000000000, 0.39818229130733, 0.38993596694726, -0.19508993437843, -0.19104961714022, 0.27189183889033, 0.26626095629775, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.51929555839097, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.00000000000000, 0.33540621099619, 0.32845997438889, -0.23868196662507, -0.23373885734945, -0.40465477419254, -0.39627436993581, 0.03194056684349, 0.24852804829932, 0.10406281520907, 0.22752003215383, 0.19164998940535, 0.39855150613730, 0.18212009644858, 0.39818229130733, 0.33540621099619, 1.00000000000000, 0.65829196968307, 0.03377788845171, 0.18973355345663, 0.03377787801250, 0.18973353027778, 0.07782687765506, 0.41673483339289, 0.06432613555760, 0.14064086581738, 0.11846790011719, 0.69577154167756, 0.17834840338088, 0.38993596694726, 0.32845997438889, 0.65829196968307, 1.00000000000000, 0.18973355345663, 0.46568275036392, 0.18973353027778, 0.46568271534951, 0.03194059360971, 0.24852817585408, 0.26133836901302, -0.11147377517723, -0.13638212535769, 0.39855159927097, 0.45736753043714, -0.19508993437843, -0.23868196662507, 0.03377788845171, 0.18973355345663, 1.00000000000000, 0.65829196968307, 0.03377753800722, 0.18973277533977, 0.07782689516847, 0.41673491190102, 0.16154551380969, -0.06890717331451, -0.08430419382570, 0.69577163288244, 0.44789543786962, -0.19104961714022, -0.23373885734945, 0.18973355345663, 0.46568275036392, 0.65829196968307, 1.00000000000000, 0.18973277533977, 0.46568157492523, 0.03194058979707, 0.24852815768492, -0.14739774471992, 0.15535813444705, -0.23121845462600, 0.39855158600482, -0.25796038297464, 0.27189183889033, -0.40465477419254, 0.03377787801250, 0.18973353027778, 0.03377753800722, 0.18973277533977, 1.00000000000000, 0.65829196968307, 0.07782689267382, 0.41673490071816, -0.09111347019718, 0.09603416103637, -0.14292698856849, 0.69577161989103, -0.25261802096777, 0.26626095629775, -0.39627436993581, 0.18973353027778, 0.46568271534951, 0.18973277533977, 0.46568157492523, 0.65829196968307, 1.00000000000000;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}