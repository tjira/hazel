#include "hf.h"

HF::ResultsRestricted HF::rscf(const System& system, Matrix D, bool print) const {
    // create all the necessary matrices, calculate ERI and initialize DIIS
    Matrix H = system.ints.T + system.ints.V, C, F; Vector eps; int nocc = system.electrons / 2;
    Tensor<4> ERI = system.ints.J - 0.5 * system.ints.J.shuffle(Array<4>{0, 3, 2, 1});
    libint2::DIIS<Matrix> diis(ropt.diis.start, ropt.diis.keep);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // calculate the Fock matrix
    if (system.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(D), Axes<2>{first, second}));
    else F = H + Integral::Coulomb(system, D);

    // calculate the energy
    double Eel = 0.5 * D.cwiseProduct(H + F).sum();

    // print the iteration header
    if (print) std::printf("\nITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= ropt.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();
        
        // calculate the Fock matrix
        if (system.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(D), Axes<2>{first, second}));
        else F = H + Integral::Coulomb(system, D);

        // exrapolate the fock matrix
        Matrix e = system.ints.S * D * F - F * D * system.ints.S;
        if (i > 1) diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, system.ints.S);

        // exteract the eigenvalues and eigenvectors and save previous values of D and E
        C = solver.eigenvectors(), eps = solver.eigenvalues();
        Matrix Dp = D; double Ep = Eel;

        // calculate the new density and energy
        D = 2 * C.leftCols(nocc) * C.leftCols(nocc).transpose();
        Eel = 0.5 * D.cwiseProduct(H + F).sum();

        // calculate the E and D errors and elapsed time
        double Eerr = std::abs(Eel - Ep), Derr = (D - Dp).norm();
        const char* elapsed = Timer::Format(Timer::Elapsed(start)).c_str();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, Eel, Eerr, Derr, elapsed, i > ropt.diis.start - 1 ? "DIIS" : "");

        // finish if covergence reached
        if (Eerr < ropt.thresh && Derr < ropt.thresh) break;
        else if (i == ropt.maxiter) {
            throw std::runtime_error("Maximum number of iterations in SCF reached.");
        }
    }

    // return the results
    return {C, D, eps, Eel + Integral::Repulsion(system), Eel, Integral::Repulsion(system)};
}
