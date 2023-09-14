#include "qdyn.h"

#define I std::complex<double>(0, 1)

Qdyn::Results Qdyn::run(System system, bool print) const {
    if (std::log2(opt.points) != (int)std::log2(opt.points)) throw std::runtime_error("NUMBER OF POINTS HAS TO BE POWER OF TWO");

    // create the output vectors
    std::vector<std::vector<CVector>> states(opt.nstates); Vector energy(opt.nstates);

    // create dx and the real and momentum space
    double dx = 2 * opt.range / (opt.points - 1);
    CVector x(opt.points), k(opt.points);
    k.fill(2 * M_PI / k.size() / dx);

    // fill the real and momentum space
    for (int i = 0; i < x.size(); i++) x(i) = 2 * i * opt.range / (opt.points - 1) - opt.range;
    for (int i = 0; i < k.size(); i++) k(i) *= i - (i < opt.points / 2 ? 0 : opt.points);

    // define the potential and the guess wfn
    CVector psi = (-(x.array() - 0.5).pow(2)).exp();
    CVector V = 0.5 * x.cwiseProduct(x);

    // create the real space and momentum space operators
    CVector K = (-0.5 * k.array().pow(2) * opt.dt).exp();
    CVector R = (-0.5 * V.array() * opt.dt).exp();

    // real time dynamics operators
    if (!opt.imaginary) {
        K = (-0.5 * I * k.array().pow(2) * opt.dt).exp();
        R = (-0.5 * I * V.array() * opt.dt).exp();
    }

    // calculate the total energy
    CVector Ek = 0.5 * Eigen::Conj(psi).array() * Eigen::fftINV(k.array().pow(2) * Eigen::fftFWD(psi).array()).array();
    CVector Ep = Eigen::Conj(psi).array() * V.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

    // loop over all states
    for (int i = 0; i < opt.nstates; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i\nITER       Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i);

        // propagate the state
        for (int j = 1; j <= opt.iters; j++) {
            // save the previous values
            CVector Dprev = psi.array().abs2(); double Eprev = E;

            // Trotter formula
            psi = R.array() * psi.array();
            psi = Eigen::fftINV(K.array() * Eigen::fftFWD(psi).array());
            psi = R.array() * psi.array();

            // subtract lower eigenstates
            for (int k = 0; k < i; k++) {
                psi = psi.array() - (Eigen::Conj(states.at(k).at(states.at(k).size() - 1)).array() * psi.array()).sum() * dx * states.at(k).at(states.at(k).size() - 1).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);

            // calculate the total energy
            Ek = 0.5 * Eigen::Conj(psi).array() * Eigen::fftINV(k.array().pow(2) * Eigen::fftFWD(psi).array()).array();
            Ep = Eigen::Conj(psi).array() * V.array() * psi.array(); E = (Ek + Ep).sum().real() * dx;

            // calculate the errors
            double Eerr = std::abs(E - Eprev), Derr = (psi.array().abs2() - Dprev.array()).abs2().sum();

            // print the iteration
            if (print) std::printf("%4d %20.14f %.2e %.2e\n", j, E, Eerr, Derr);

            // append the wave function
            states.at(i).push_back(psi);

            // end the loop if converged
            if (Eerr < opt.thresh && Derr < opt.thresh) break;
            else if (j == opt.iters && opt.imaginary) {
                throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN ITP REACHED.");
            }
        }

        // assign the energy
        energy(i) = E;
    }

    // return the results
    return {states, energy, x};
}
