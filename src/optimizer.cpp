#include "optimizer.h"

template <class M>
Optimizer<M>::Results Optimizer<M>::optimize(const System& system, bool print) const {
    // run the HF optimizer
    if constexpr (std::is_same_v<HF, M>) {
        auto egfunc = [](System& system, Data data) {
            system.dints.dS = Integral::dOverlap(system), system.dints.dT = Integral::dKinetic(system);
            system.dints.dV = Integral::dNuclear(system), system.dints.dJ = Integral::dCoulomb(system);
            system.ints.S = Integral::Overlap(system), system.ints.T = Integral::Kinetic(system);
            system.ints.V = Integral::Nuclear(system), system.ints.J = Integral::Coulomb(system);
            return Gradient<HF>(HF(data).rscf(system, false)).get(system, false);
        };
        return optimize(system, egfunc, print);

    // run the MP optimizer
    } else if constexpr (std::is_same_v<MP, M>) {
        auto egfunc = [](System& system, Data data) {
            system.ints.J = Integral::Coulomb(system);
            system.ints.S = Integral::Overlap(system);
            system.ints.T = Integral::Kinetic(system);
            system.ints.V = Integral::Nuclear(system);
            data = HF(data).rscf(system, false);
            HF::ResultsRestricted rhfres = {data.hf.C, data.hf.D, data.hf.eps, data.hf.E, data.hf.E - Integral::Repulsion(system), Integral::Repulsion(system)};
            Tensor<4> Jmo = Transform::Coulomb(system.ints.J, data.hf.C);
            data.mp.Ecorr = MP({rhfres}).mp2(system, Jmo, false).Ecorr;
            return Gradient<MP>(data).get(system, false);
        };
        return optimize(system, egfunc, print);
    }
}

template <class M>
Optimizer<M>::Results Optimizer<M>::optimize(const System& system, const std::function<Data(System&, Data)>& egfunc, bool print) const {
    System tempsys = system;
    // create the output and perform the SCF and gradient calculation
    Data output = data; output = egfunc(tempsys, data); Matrix G; double thresh, E;

    // assign the correct threshond
    if constexpr (std::is_same_v<MP, M>) G = output.mp.grad.G, thresh = data.mp.opt.thresh, E = output.hf.E + output.mp.Ecorr;
    if constexpr (std::is_same_v<HF, M>) G = output.hf.grad.G, thresh = data.hf.opt.thresh, E = output.hf.E;

    // print the header
    if (print) std::printf("ITER        E [Eh]         |GRAD|      TIME\n");

    // print the initial state info
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");


    // move the system while gradient is big
    for (int i = 1; G.norm() > thresh; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // move the system
        tempsys.move(-G);

        // perform HF and calculate gradient
        output = egfunc(tempsys, output);

        // extract the correct energy and gradient
        if constexpr (std::is_same_v<MP, M>) G = output.mp.grad.G, E = output.hf.E + output.mp.Ecorr;
        if constexpr (std::is_same_v<HF, M>) G = output.hf.grad.G, E = output.hf.E;

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return results
    return {tempsys, G};
}

template class Optimizer<HF>;
template class Optimizer<MP>;
