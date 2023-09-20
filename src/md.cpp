#include "md.h"

MD::Results MD::run(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print) const {
    // calculate the initial energy and gradient and print the header with the initial state info
    auto[E, G] = egfunc(system); if (print) std::printf(" ITER         E [Eh]         |GRAD|      TIME\n");
    if (print) std::printf("%6d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");

    // create position, velocity and mass matrices
    Matrix v(system.atoms.size(), 3), a(system.atoms.size(), 3), m(system.atoms.size(), 3);

    // fill the mass matrix
    for(size_t i = 0; i < system.atoms.size(); i++) {
        m.row(i) = [](double m) {return Vector::Constant(3, m);}(masses.at(system.atoms.at(i).atomic_number));
    }

    // write the initial geometry
    system.save(opt.output);

    // move the system while gradient is big
    for (int i = 0; i < opt.iters; i++) {
        // start the timer and store the previous v and a
        auto start = Timer::Now();
        Matrix vp = v, ap = a;

        // calculate the velocity and accceleration
        a = -G.array() / m.array(); v = vp + (ap + a) * opt.step / 2;

        // move the system and calculate the next energy with gradient
        system.move(opt.step * (v + 0.5 * a * opt.step)), std::tie(E, G) = egfunc(system);

        // write the current geometry
        system.save(opt.output, std::ios::app);

        // print the iteration info
        if (print) std::printf("%6d %20.14f %.2e %s\n", i + 1, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return the results
    return {};
}
