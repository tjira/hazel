#include "qdyn.h"

Qdyn::Results Qdyn::run(System system, bool print) const {
    if (std::log2(opt.points) != (int)std::log2(opt.points)) throw std::runtime_error("NUMBER OF POINTS HAS TO BE POWER OF TWO");

    // create the output vectors
    std::vector<CVector> states;
    Vector energy(opt.nstates);

    // create dx and the real and momentum space
    double dx = 2 * opt.range / (opt.points - 1);
    CVector x(opt.points), k(opt.points);
    k.fill(2 * M_PI / k.size() / dx);

    // fill the real and momentum space
    for (int i = 0; i < x.size(); i++) x(i) = 2 * i * opt.range / (opt.points - 1) - opt.range;
    for (int i = 0; i < k.size(); i++) k(i) *= i - (i < opt.points / 2 ? 0 : opt.points);

    // define the potential
    CVector V = 0.5 * x.cwiseProduct(x);

    // create the real space and momentum space operators
    CVector K = (-0.5 * k.array().pow(2) * opt.dt).exp();
    CVector R = (-0.5 * V.array() * opt.dt).exp();

    // loop over all states
    for (int i = 0; i < opt.nstates; i++) {
        // create the guess wave function
        CVector psi = (-(x.array() - 0.5).pow(2)).exp();

        // propagate the wavefunction
        for (int j = 0; j < opt.iters; j++) {
            // Trotter formula
            psi = R.array() * psi.array();
            psi = Eigen::fftINV(K.array() * Eigen::fftFWD(psi).array());
            psi = R.array() * psi.array();

            // subtract lower eigenstates
            for (int k = 0; k < i; k++) {
                psi = psi.array() - (Eigen::Conj(states.at(k)).array() * psi.array()).sum() * dx * states.at(k).array();
            }

            // normalize the wavefunction
            psi = psi.array() / std::sqrt(psi.array().abs2().sum() * dx);
        }

        // append the wave function
        states.push_back(psi);

        // calculate the total energy
        CVector Ek = 0.5 * Eigen::Conj(psi).array() * Eigen::fftINV(k.array().pow(2) * Eigen::fftFWD(psi).array()).array();
        CVector Ep = Eigen::Conj(psi).array() * V.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

        // assign the energy
        energy(i) = E;
    }

    // return the results
    return {states, energy, x};
}
