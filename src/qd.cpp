#include "qd.h"

#define I std::complex<double>(0, 1)

QD::Results QD::run(System, bool print) const {
    // create the output energy vector
    Vector energy(opt.nstates);

    // calculate the potential lines and create all the vectors
    std::ifstream stream(opt.potfile); std::string line; int lines; std::getline(stream, line);
    for (lines = 0; std::getline(stream, line); lines++){}; stream.clear(), stream.seekg(0); 
    CVector x(lines), k(lines), V(lines); std::getline(stream, line);

    // read the potential
    for (int i = 0; i < lines; i++) {
        std::getline(stream, line); std::istringstream iss(line); iss >> x(i); iss >> V(i);
    }

    // calculate dx and create the momentum space and guess wfn
    double dx = (x(1) - x(0)).real(); k.fill(2 * M_PI / k.size() / dx);
    CVector psi = (-(x.array() - 0.8).pow(2)).exp(); stream.close();
    std::vector<std::vector<CVector>> states(opt.nstates, {psi}); 

    // fill the real and momentum space
    for (int i = 0; i < k.size(); i++) k(i) *= i - (i < x.size() / 2 ? 0 : x.size());

    // create the momentum space operators
    CVector K = (-0.5 * k.array().pow(2) * opt.dt).exp();
    CVector R = (-0.5 * V.array() * opt.dt).exp();

    // real time dynamics operators
    if (!opt.imaginary) {
        K = (-0.5 * I * k.array().pow(2) * opt.dt).exp();
        R = (-0.5 * I * V.array() * opt.dt).exp();

        // create imaginary options
        Options imopt = opt; imopt.imaginary = true;
        imopt.iters = 1000;

        // optimize the wavefunctions
        states = QD(imopt).run(System(), false).states;

        // delete optimization steps
        for (size_t i = 0; i < states.size(); i++) {
            states.at(i) = {states.at(i).at(states.at(i).size() - 1)};
        }
    }

    // calculate the total energy
    CVector Ek = 0.5 * Eigen::Conj(psi).array() * Eigen::fftINV(k.array().pow(2) * Eigen::fftFWD(psi).array()).array();
    CVector Ep = Eigen::Conj(psi).array() * V.array() * psi.array(); double E = (Ek + Ep).sum().real() * dx;

    // loop over all states
    for (int i = 0; i < opt.nstates; i++) {
        // print the iteration header
        if (print) std::printf("%sSTATE %i\nITER       Eel [Eh]         |dE|     |dD|\n", i ? "\n" : "", i);

        // assign the psi WFN to the correct state
        psi = states.at(i).at(states.at(i).size() - 1);

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
            if (opt.imaginary && Eerr < opt.thresh && Derr < opt.thresh) break;
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
