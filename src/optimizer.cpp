#include "optimizer.h"

template <class M>
Data Optimizer<M>::optimize(bool print) const {
    // run the HF optimizer
    if constexpr (std::is_same_v<HF, M>) {
        auto egfunc = [](Data data) {
            data.ints.dS = Integral::dOverlap(data.system), data.ints.dT = Integral::dKinetic(data.system);
            data.ints.dV = Integral::dNuclear(data.system), data.ints.dJ = Integral::dCoulomb(data.system);
            data.ints.S = Integral::Overlap(data.system), data.ints.T = Integral::Kinetic(data.system);
            data.ints.V = Integral::Nuclear(data.system), data.ints.J = Integral::Coulomb(data.system);
            return Gradient<HF>(HF(data).scf(false)).get(false);
        };
        return optimize(egfunc, print);

    // run the MP optimizer
    } else if constexpr (std::is_same_v<MP, M>) {
        auto egfunc = [](Data data) {
            data.ints.S = Integral::Overlap(data.system), data.ints.T = Integral::Kinetic(data.system);
            data.ints.V = Integral::Nuclear(data.system), data.ints.J = Integral::Coulomb(data.system);
            data = HF(data).scf(false); data.intsmo.J = Transform::Coulomb(data.ints.J, data.hf.C);
            return Gradient<MP>(MP(data).mp2(false)).get(false);
        };
        return optimize(egfunc, print);
    }
}

template <class M>
Data Optimizer<M>::optimize(const std::function<Data(Data)>& egfunc, bool print) const {
    // create the output and perform the SCF and gradient calculation
    Data output = data; output = egfunc(data); Matrix G; double thresh, E;

    // assign the correct threshond
    if constexpr (std::is_same_v<MP, M>) G = output.mp.grad.G, thresh = data.mp.opt.thresh, E = output.hf.E + output.mp.Ecorr;
    if constexpr (std::is_same_v<HF, M>) G = output.hf.grad.G, thresh = data.hf.opt.thresh, E = output.hf.E;

    // print the header
    if (print) std::printf("ITER        E [Eh]         |GRAD|      TIME\n");

    // print the initial state info
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");

    // move the data.system while gradient is big
    for (int i = 1; G.norm() > thresh; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // move the system
        output.system.move(-G);

        // perform HF and calculate gradient
        output = egfunc(output);

        // extract the correct energy and gradient
        if constexpr (std::is_same_v<MP, M>) G = output.mp.grad.G, E = output.hf.E + output.mp.Ecorr;
        if constexpr (std::is_same_v<HF, M>) G = output.hf.grad.G, E = output.hf.E;

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return results
    return output;
}

template class Optimizer<HF>;
template class Optimizer<MP>;
