#include "roothaan.h"

Roothaan::Roothaan(const System& system, int maxiter, double thresh, const std::pair<int, int>& diis) : diis(diis.first - 2, diis.second), thresh(thresh), maxiter(maxiter), system(system) {}

Matrix Roothaan::gradient(const Integrals& ints, const Matrix& C, const Vector& eps) const {
    // extract the useful stuff from the calculated integrals and define all the contractio axes
    Tensor<3> dS1 = ints.dS.slice<Index<3>, Index<3>>({0, 0, 0}, {ints.dS.dimension(0), ints.dS.dimension(1), 3});
    Tensor<3> dT1 = ints.dT.slice<Index<3>, Index<3>>({0, 0, 0}, {ints.dT.dimension(0), ints.dT.dimension(1), 3});
    Tensor<3> dV1 = ints.dV.slice<Index<3>, Index<3>>({0, 0, 0}, {ints.dV.dimension(0), ints.dV.dimension(1), 3});
    Pair first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = system.electrons / 2;

    // define the density, weighed density and gradient matrix
    Matrix G = Matrix::Zero(system.atoms.size(), 3), D = 2 * C.leftCols(nocc) * C.leftCols(nocc).transpose();
    Tensor<2> W(C.rows(), C.cols()); auto atom2shell = system.shells.atom2shell(system.atoms);

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * C.leftCols(nocc).row(i).cwiseProduct(C.leftCols(nocc).row(j)) * eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI tensor
    Tensor<3> dERI = (ints.dJ - 0.5 * ints.dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(D), Axes<2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.shells.at(shell).size();

        // define the Hcore derivative and atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = ints.dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {D.rows(), D.cols(), 3});
        Eigen::array<Eigen::Index, 3> Soff = {si, 0, 0}, Sext = {ss, D.cols(), 3};
        Eigen::array<Eigen::Index, 2> Doff = {si, 0}, Dext = {ss, D.cols()};

        // fill the Hcore derivative
        dHcore.slice<Index<3>, Index<3>>({0, si, 0}, {D.rows(), ss, 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext).shuffle(Index<3>{1, 0, 2});
        dHcore.slice<Index<3>, Index<3>>({si, 0, 0}, {ss, D.cols(), 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext);

        // contract the tensors and add them to the gradient
        G.row(i) += 2 * toVector(dERI.slice(Soff, Sext).contract(toTensor(D).slice(Doff, Dext), Axes<2>{third, fourth}));
        G.row(i) -= 2 * toVector(dS1.slice(Soff, Sext).contract(W.slice(Doff, Dext), Axes<2>{third, fourth}));
        G.row(i) += toVector(dHcore.contract(toTensor(D), Axes<2>{third, fourth}));
    }

    // return the gradient
    return G + Integral::dRepulsion(system);
}

std::tuple<System, Integrals, Matrix, Matrix> Roothaan::optimize(Integrals ints) const {
    // define Densty matrix, perform HF calculation and calculate the analytical gradient
    Matrix D = Matrix::Zero(ints.S.rows(), ints.S.cols()); auto[C, eps, E] = scf(ints, D, false);
    Matrix G = gradient(ints, C, eps); System optsys = system;

    // print the header
    std::printf("ITER       E [Eh]         |dG|        TIME    \n");

    // print the initial state info
    std::printf("%4d %20.14f %.2e %s\n", 0, E, G.norm(), "00:00:00.000");

    // move the system while gradient is big
    for (int i = 1; G.norm() > thresh; i++) {
        // start the timer and move the system
        Timer::Timepoint start = Timer::Now();
        optsys.move(-G);

        // calculate integrals for HF method
        ints.S = Integral::Overlap(optsys);
        ints.T = Integral::Kinetic(optsys);
        ints.V = Integral::Nuclear(optsys);
        ints.J = Integral::Coulomb(optsys);

        // calculate integral derivatives for the gradient
        ints.dS = Integral::dOverlap(optsys);
        ints.dT = Integral::dKinetic(optsys);
        ints.dV = Integral::dNuclear(optsys);
        ints.dJ = Integral::dCoulomb(optsys);

        // perform HF and calculate gradient
        Roothaan rooth(optsys, maxiter, thresh, {diis.start, diis.keep});
        std::tie(C, eps, E) = rooth.scf(ints, D, false);
        G = rooth.gradient(ints, C, eps);

        // print the iteration info
        std::printf("%4d %20.14f %.2e %s\n", i, E, G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return stuff
    return {optsys, ints, D, G};
}

std::tuple<Matrix, Vector, double> Roothaan::scf(const Integrals& ints, Matrix& D, bool print) const {
    // create all the necessary matrices, calculate ERI and initialize DIIS
    Tensor<4> ERI = ints.J - 0.5 * ints.J.shuffle(Array<4>{0, 3, 2, 1});
    libint2::DIIS<Matrix> diis(this->diis.start, this->diis.keep);
    Matrix eps(D.rows(), 1), C(D.rows(), D.cols());
    Matrix H = ints.T + ints.V; double E = 0;
    int nocc = system.electrons / 2;

    // print the iteration header
    if (print) std::printf("ITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= maxiter; i++) {
        // define the contraction axes and start the timer
        Eigen::IndexPair<int> first(2, 0), second(3, 1);
        Timer::Timepoint start = Timer::Now();

        // calculate the fock matrix and extrapolate it
        Matrix F = H + toMatrix(ERI.contract(toTensor(D), Axes<2>{first, second}));
        Matrix e = ints.S * D * F - F * D * ints.S;
        if (i > 1) diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, ints.S);
        C = solver.eigenvectors(), eps = solver.eigenvalues();
        Matrix Dp = D; double Ep = E;

        // calculate the new density and energy
        D = 2 * C.leftCols(nocc) * C.leftCols(nocc).transpose();
        E = 0.5 * D.cwiseProduct(H + F).sum();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, E, std::abs(E - Ep), (D - Dp).norm(), Timer::Format(Timer::Elapsed(start)).c_str(), i > this->diis.start + 1 ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(E - Ep) < thresh && (D - Dp).norm() < thresh) break;
    }

    // return the results
    return {C, eps, E + Integral::Repulsion(system)};
}
