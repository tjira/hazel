#include "md.h"

template <typename F>
void MD::run(System system, const F& egfunc, bool print) const {
    // define energy and gradients
    Vector E(1); std::vector<Matrix> G;

    // calculate the initial energy and gradient
    if constexpr (std::is_same_v<F, std::function<std::tuple<Vector, std::vector<Matrix>>(System&, const std::vector<int>&)>>) {
        std::tie(E, G) = egfunc(system, {opt.state});
    } else if constexpr (std::is_same_v<F, std::function<std::tuple<double, Matrix>(System&)>>) {
        auto[E0, G0] = egfunc(system); E << E0, G.push_back(G0);
    }

    // print the header
    if (print) std::printf("\n ITER  TIME [fs] S");
    if (print) for (int i = 0; i < E.size(); i++) std::printf("        E [Eh]       ");
    if (print) std::printf("       KIN [Eh]              TK [K]         |GRAD|      TIME\n");

    // get the degrees of freedom
    double Nf = system.atoms.size() * 3 - 6;

    // create position, velocity and mass matrices
    Matrix v(system.atoms.size(), 3), a(system.atoms.size(), 3), m(system.atoms.size(), 3);

    // fill the mass matrix
    for(size_t i = 0; i < system.atoms.size(); i++) {
        m.row(i) = [](double m) {return Vector::Constant(3, m);}(masses.at(system.atoms.at(i).atomic_number));
    }

    // write the initial geometry and define state
    system.save(opt.output); int state = opt.state;

    // print the zeroth iteration
    if (print) std::printf("%6d %9.4f %d", 0, 0.0, opt.state);
    if (print) for (int i = 0; i < E.size(); i++) std::printf(" %20.14f", E(i));
    if (print) std::printf(" %20.14f %20.14f %.2e %s\n", 0.0, 0.0, G.at(0).norm(), "00:00:00.000");

    // move the system while gradient is big
    for (int i = 0; i < opt.iters; i++) {
        // start the timer and store the previous v and a
        auto start = Timer::Now();
        Matrix vp = v, ap = a;

        // calculate the velocity and accceleration
        a = -G.at(0).array() / m.array(); v = vp + (ap + a) * opt.step / 2;

        // calculate the kinetic energy and temperature before thermostatting
        double Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); double T = 2 * Ekin / Nf;

        // apply the berendsen thermostat
        if (opt.berendsen.temp > 0) v.array() *= std::sqrt(1 + (opt.berendsen.temp / T - 1) * opt.step / opt.berendsen.tau);

        // calculate the kinetic energy and temperature after thermostatting
        Ekin = 0.5 * (m.array() * v.array() * v.array()).sum(); T = 2 * Ekin / Nf;

        // move the system
        system.move(opt.step * (v + 0.5 * a * opt.step));

        // write the current geometry
        system.save(opt.output, std::ios::app);

        // calculate the next energy with gradient
        if constexpr (std::is_same_v<F, std::function<std::tuple<Vector, std::vector<Matrix>>(System&, const std::vector<int>&)>>) {
            std::tie(E, G) = egfunc(system, {opt.state});
        } else if constexpr (std::is_same_v<F, std::function<std::tuple<double, Matrix>(System&)>>) {
            auto[E0, G0] = egfunc(system); E << E0, G = {G0};
        }

        // print the iteration info
        if (print) std::printf("%6d %9.4f %d", i + 1, AU2FS * opt.step * (i + 1), state);
        if (print) for (int i = 0; i < E.size(); i++) std::printf(" %20.14f", E(i));
        if (print) std::printf(" %20.14f %20.14f %.2e %s\n", Ekin, T / BOLTZMANN, G.at(0).norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }
}

template void MD::run(System, const std::function<std::tuple<Vector, std::vector<Matrix>>(System&, const std::vector<int>&)>&, bool) const;
template void MD::run(System, const std::function<std::tuple<double, Matrix>(System&)>&, bool) const;
