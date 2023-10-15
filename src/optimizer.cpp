#include "optimizer.h"

System Optimizer::optimize(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print) const {
    // calculate the initial energy and gradient and print the header with the initial state info
    auto[E, G] = egfunc(system); if (print) std::printf("\nITER        E [Eh]         |GRAD|      TIME\n");
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");

    // move the system while gradient is big
    for (int i = 1; G.norm() > opt.thresh; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // move the system and calculate energy with gradient
        system.move(-G), std::tie(E, G) = egfunc(system);

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return the system
    return system;
}
