#include "qdyn.h"
#include "eigen.h"

Qdyn::Results Qdyn::run(System system, bool print) const {
    if (std::log2(opt.points) != (int)std::log2(opt.points)) throw std::runtime_error("NUMBER OF POINTS HAS TO BE POWER OF TWO");

    // create the output vectors
    std::vector<CVector> states;
    Vector energy(opt.nstates);

    // create dx and the real and momentum space
    double dx = 2 * opt.range / (opt.points - 1);
    CVector x(opt.points), k(opt.points);

    // fil the real space
    for (int i = 0; i < x.size(); i++) {
        x(i) = -opt.range + i * dx;
    }

    // define the potential
    CVector V = 0.5 * x.cwiseProduct(x);

    // fill the momentum space
    for (int i = 0; i < k.size(); i++) {
        if (int half = k.size() / 2; i < half) k(i) = i;
        else k(i) = -half + i - half;
    }
    k /= dx * k.size() / M_PI / 2;

    // create the real space and momentum space operators
    CVector R = (-0.5 * V.array() * opt.dt).exp(), K = (-0.5 * k.array().pow(2) * opt.dt).exp();

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

    // OUTPUT
    std::ofstream file("wavefunction.dat");
    file << "# x";
    for (int i = 0; i < opt.nstates; i++) {
        file << " state" << i << ".real " << "state" << i << ".imag";
    }
    file << "\n";
    for (int j = 0; j < x.size(); j++) {
        file << x(j).real();
        for (int i = 0; i < opt.nstates; i++) {
            file << " " << states.at(i)(j).real() << " " << states.at(i)(j).imag();
        }
        file << "\n";
    }

    return {states, energy, x};
}
