#include "hf.h"

HF::ResultsRestricted HF::rscf(const System& system, Matrix D, bool print) const {
    // calculate integrals if not provided
    if (!system.ints.J.size() && !ropt.nocoulomb) const_cast<System&>(system).ints.J = Integral::Coulomb(system);
    if (!system.ints.S.size()) const_cast<System&>(system).ints.S = Integral::Overlap(system);
    if (!system.ints.T.size()) const_cast<System&>(system).ints.T = Integral::Kinetic(system);
    if (!system.ints.V.size()) const_cast<System&>(system).ints.V = Integral::Nuclear(system);

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
        std::string elapsed = Timer::Format(Timer::Elapsed(start));
        double Eerr = std::abs(Eel - Ep), Derr = (D - Dp).norm();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, Eel, Eerr, Derr, elapsed.c_str(), i > ropt.diis.start - 1 ? "DIIS" : "");

        // finish if covergence reached
        if (Eerr < ropt.thresh && Derr < ropt.thresh) break;
        else if (i == ropt.maxiter) {
            throw std::runtime_error("MAXIMUM NUMBER OF ITERATIONS IN SCF REACHED.");
        }
    }

    // return the results
    return {C, D, eps, Eel + Integral::Repulsion(system), Eel, Integral::Repulsion(system)};
}

HF::ResultsUnrestricted HF::uscf(const System& system, Matrix D, bool print) const {
    // calculate integrals if not provided
    if (!system.ints.J.size() && !uopt.nocoulomb) const_cast<System&>(system).ints.J = Integral::Coulomb(system);
    if (!system.ints.S.size()) const_cast<System&>(system).ints.S = Integral::Overlap(system);
    if (!system.ints.T.size()) const_cast<System&>(system).ints.T = Integral::Kinetic(system);
    if (!system.ints.V.size()) const_cast<System&>(system).ints.V = Integral::Nuclear(system);

    // create all the necessary matrices, calculate K
    Tensor<4> K = system.ints.J.shuffle(Array<4>{0, 3, 2, 1}); Matrix Da = D, Db = D;
    Matrix H = system.ints.T + system.ints.V, Ca, Cb, Fa, Fb; Vector epsa, epsb;

    // initialize DIIS and contraction axes
    libint2::DIIS<Matrix> diisa(uopt.diis.start, uopt.diis.keep);
    libint2::DIIS<Matrix> diisb(uopt.diis.start, uopt.diis.keep);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // calculate number of alpha and beta electrons
    int na = (system.electrons + system.multi - 1) / 2;
    int nb = (system.electrons - system.multi + 1) / 2;

    // create coulomb and exchange matrices for both spins
    Matrix Ja = toMatrix(system.ints.J.contract(toTensor(Da), Axes<2>{first, second}));
    Matrix Jb = toMatrix(system.ints.J.contract(toTensor(Db), Axes<2>{first, second}));
    Matrix Ka = toMatrix(K.contract(toTensor(Da), Axes<2>{first, second}));
    Matrix Kb = toMatrix(K.contract(toTensor(Db), Axes<2>{first, second}));

    // calculate the Fock matrix for alpha electrons
    if (system.ints.J.size()) Fa = H + 0.5 * (Ja + Jb - Ka);
    else throw std::runtime_error("COULOMB TENSOR HAS TO BE CALCULATED FOR UHF");

    // calculate the Fock matrix for beta electrons
    if (system.ints.J.size()) Fb = H + 0.5 * (Ja + Jb - Kb);
    else throw std::runtime_error("COULOMB TENSOR HAS TO BE CALCULATED FOR UHF");

    // calculate the electron energy
    double Eel = 0.25 * (Da.cwiseProduct(H + Fa).sum() + Db.cwiseProduct(H + Fb).sum());

    // print the iteration header
    if (print) std::printf("\nITER       Eel [Eh]         |dE|    |dD_a|   |dD_b|      TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= uopt.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // create coulomb and exchange matrices for both spins
        Ja = toMatrix(system.ints.J.contract(toTensor(Da), Axes<2>{first, second}));
        Jb = toMatrix(system.ints.J.contract(toTensor(Db), Axes<2>{first, second}));
        Ka = toMatrix(K.contract(toTensor(Da), Axes<2>{first, second}));
        Kb = toMatrix(K.contract(toTensor(Db), Axes<2>{first, second}));
        
        // calculate the Fock matrix for alpha electrons
        if (system.ints.J.size()) Fa = H + 0.5 * (Ja + Jb - Ka);
        else throw std::runtime_error("COULOMB TENSOR HAS TO BE CALCULATED FOR UHF");

        // calculate the Fock matrix for beta electrons
        if (system.ints.J.size()) Fb = H + 0.5 * (Ja + Jb - Kb);
        else throw std::runtime_error("COULOMB TENSOR HAS TO BE CALCULATED FOR UHF");

        // exrapolate the fock matrix
        Matrix e = system.ints.S * Da * Fa - Fa * Da * system.ints.S + system.ints.S * Db * Fb - Fb * Db * system.ints.S;
        if (i > 1) diisa.extrapolate(Fa, e), diisb.extrapolate(Fb, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solvera(Fa, system.ints.S);
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solverb(Fb, system.ints.S);

        // exteract the eigenvalues and eigenvectors and save previous values of D and E
        Ca = solvera.eigenvectors(), epsa = solvera.eigenvalues();
        Cb = solverb.eigenvectors(), epsb = solverb.eigenvalues();
        Matrix Dap = Da, Dbp = Db; double Ep = Eel;

        // calculate the density matrices
        Da = 2 * Ca.leftCols(na) * Ca.leftCols(na).transpose();
        Db = 2 * Cb.leftCols(nb) * Cb.leftCols(nb).transpose();

        // calculate the new energy
        Eel = 0.25 * (Da.cwiseProduct(H + Fa).sum() + Db.cwiseProduct(H + Fb).sum());

        // calculate the E and D errors and elapsed time
        double Eerr = std::abs(Eel - Ep), Daerr = (Da - Dap).norm(), Dberr = (Db - Dbp).norm();
        std::string elapsed = Timer::Format(Timer::Elapsed(start));

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %.2e %s %s\n", i, Eel, Eerr, Daerr, Dberr, elapsed.c_str(), i > uopt.diis.start - 1 ? "DIIS" : "");

        // finish if covergence reached
        if (Eerr < uopt.thresh && Daerr < uopt.thresh) break;
        else if (i == uopt.maxiter) {
            throw std::runtime_error("Maximum number of iterations in SCF reached.");
        }
    }

    // return the results
    return {Ca, Cb, Da, Db, Eel + Integral::Repulsion(system), Eel, Integral::Repulsion(system), epsa, epsb};
}
