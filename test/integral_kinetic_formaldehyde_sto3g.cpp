#include "../include/integral.h"

int test_integral_kinetic_formaldehyde_sto3g(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "STO-3G", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix I = Integral::Kinetic(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Iexp(system.shells.nbf(), system.shells.nbf()); Iexp << 29.00320406467810, -0.16801096113783, -0.00000000000000, 0.00000000000000, -0.00000000000000, -0.00005270259706, -0.00520720038093, 0.00303679945472, 0.00007682799404, 0.00034965860204, -0.00166755594240, -0.00166755691351, -0.16801096113783, 0.80812790277365, -0.00000000000000, -0.00000000000000, -0.00000000000000, -0.01914510978177, 0.06047223668844, -0.21778004565418, -0.00550961770730, -0.02507530294685, -0.01730891190031, -0.01730891886778, -0.00000000000000, -0.00000000000000, 2.52873122631647, 0.00000000000000, 0.00000000000000, -0.02227510635432, 0.16312111608817, -0.38841156785571, -0.01226039274580, -0.05579934555915, -0.00600095205968, -0.00670627521602, 0.00000000000000, -0.00000000000000, 0.00000000000000, 2.52873122631647, 0.00000000000000, -0.00056353794965, 0.00412680136481, -0.01226039274580, 0.09589788593526, -0.00141166772844, 0.00044090892879, -0.00076238894864, -0.00000000000000, -0.00000000000000, 0.00000000000000, 0.00000000000000, 2.52873122631647, -0.00256476684230, 0.01878184656750, -0.05579934555915, -0.00141166772844, 0.08978329718513, -0.00392668161892, 0.00246356407247, -0.00005270259706, -0.01914510978177, -0.02227510635432, -0.00056353794965, -0.00256476684230, 15.89112181239595, -0.08588999412233, -0.00000000000000, -0.00000000000000, -0.00000000000000, -0.01089842593288, -0.01089842619719, -0.00520720038093, 0.06047223668844, 0.16312111608817, 0.00412680136481, 0.01878184656750, -0.08588999412233, 0.47224999256911, 0.00000000000000, 0.00000000000000, 0.00000000000000, 0.10318534333726, 0.10318534056587, 0.00303679945472, -0.21778004565418, -0.38841156785571, -0.01226039274580, -0.05579934555915, -0.00000000000000, 0.00000000000000, 1.47772809059759, 0.00000000000000, 0.00000000000000, 0.11347266065065, 0.15959127270477, 0.00007682799404, -0.00550961770730, -0.01226039274580, 0.09589788593526, -0.00141166772844, -0.00000000000000, 0.00000000000000, 0.00000000000000, 1.47772809059759, 0.00000000000000, -0.03588520840724, 0.04279348445686, 0.00034965860204, -0.02507530294685, -0.05579934555915, -0.00141166772844, 0.08978329718513, -0.00000000000000, 0.00000000000000, 0.00000000000000, 0.00000000000000, 1.47772809059759, 0.22463615283141, -0.19319558676006, -0.00166755594240, -0.01730891190031, -0.00600095205968, 0.00044090892879, -0.00392668161892, -0.01089842593288, 0.10318534333726, 0.11347266065065, -0.03588520840724, 0.22463615283141, 0.76003187992239, -0.00755817890009, -0.00166755691351, -0.01730891886778, -0.00670627521602, -0.00076238894864, 0.00246356407247, -0.01089842619719, 0.10318534056587, 0.15959127270477, 0.04279348445686, -0.19319558676006, -0.00755817890009, 0.76003187992239;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL: " << I << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED INTEGRAL NORM: " << I.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL: " << Iexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED INTEGRAL NORM: " << Iexp.norm() << std::endl;

    // return success or failure based on the error
    return (I - Iexp).norm() > 1e-8;
}