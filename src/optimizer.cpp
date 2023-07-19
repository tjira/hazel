#include "optimizer.h"

template <class M>
System Optimizer<M>::optimize(const System& system, bool print) const {
    // run the HF optimizer
    if constexpr (std::is_same_v<HF, M>) {
        auto egfunc = [this](System& system) {
            system.dints.dS = Integral::dOverlap(system), system.dints.dT = Integral::dKinetic(system);
            system.dints.dV = Integral::dNuclear(system), system.dints.dJ = Integral::dCoulomb(system);
            system.ints.J = Integral::Coulomb(system), system.ints.S = Integral::Overlap(system);
            system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
            HF::ResultsRestricted rhfres = HF(ropt.rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false);
            Gradient<HF>::OptionsRestricted gradrhfopt = {ropt.rhfopt, rhfres, ropt.gradopt.numerical, ropt.gradopt.step};
            return std::tuple{rhfres.E, Gradient<HF>(gradrhfopt).get(system, false)};
        };
        return optimize(system, egfunc, print);

    // run the MP optimizer
    } else if constexpr (std::is_same_v<MP, M>) {
        auto egfunc = [this](System& system) {
            system.ints.J = Integral::Coulomb(system), system.ints.S = Integral::Overlap(system);
            system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
            HF::ResultsRestricted rhfres = HF(ropt.rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false);
            Gradient<MP>::OptionsRestricted gradrmpopt = {ropt.rhfopt, rhfres, ropt.gradopt.numerical, ropt.gradopt.step};
            Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);
            double Ecorr = MP({rhfres}).mp2(system, Jmo, false).Ecorr;
            return std::tuple{rhfres.E + Ecorr, Gradient<MP>(gradrmpopt).get(system, false)};
        };
        return optimize(system, egfunc, print);
    }
}

template <class M>
System Optimizer<M>::optimize(const System& system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print) const {
    System tempsys = system;
    // create the output and perform the SCF and gradient calculation
    auto[E, G] = egfunc(tempsys);

    // print the header
    if (print) std::printf("ITER        E [Eh]         |GRAD|      TIME\n");

    // print the initial state info
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");


    // move the system while gradient is big
    for (int i = 1; G.norm() > ropt.thresh; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // move the system
        tempsys.move(-G);

        // perform HF and calculate gradient
        std::tie(E, G) = egfunc(tempsys);

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return results
    return tempsys;
}

template class Optimizer<HF>;
template class Optimizer<MP>;
