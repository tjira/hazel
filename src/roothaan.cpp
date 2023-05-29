#include "roothaan.h"

Roothaan::Roothaan(const System& system, int maxiter, double thresh, const std::pair<int, int>& diis) : diis(diis.first - 2, diis.second), thresh(thresh), maxiter(maxiter), system(system) {}

Matrix Roothaan::gradient(const Tensor<3>& dT, const Tensor<3>& dV, const Tensor<5>& dJ, const Tensor<3>& dS, const Matrix& C, const Vector& eps) const {
    // extract map between atoms and shells and calculate the density matrix
    Matrix D = 2 * C.leftCols(system.nocc) * C.leftCols(system.nocc).transpose();
    auto atom2shell = system.shells.atom2shell(system.atoms);

    // extract the useful stuff from the calculated integrals and define all the contractio axes
    Tensor<3> dS1 = dS.slice<Index<3>, Index<3>>({0, 0, 0}, {D.rows(), D.cols(), 3});
    Tensor<3> dT1 = dT.slice<Index<3>, Index<3>>({0, 0, 0}, {D.rows(), D.cols(), 3});
    Tensor<3> dV1 = dV.slice<Index<3>, Index<3>>({0, 0, 0}, {D.rows(), D.cols(), 3});
    Pair first(2, 0), second(3, 1), third(0, 0), fourth(1, 1);

    // calculate the derivative of the ERI tensor and define the weighed density and gradient matrix
    Tensor<3> dERI = (dJ - 0.5 * dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(D), Axes<2>{first, second});
    Tensor<2> W(D.rows(), D.cols()); Matrix G = Matrix::Zero(system.atoms.size(), 3);

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * C.leftCols(system.nocc).row(i).cwiseProduct(C.leftCols(system.nocc).row(j)) * eps.topRows(system.nocc);
        }
    }

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.shells.at(shell).size();

        // define the Hcore derivative and atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {D.rows(), D.cols(), 3});
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
    return G;
}

std::tuple<Matrix, Vector, double> Roothaan::scf(const Matrix& H, const Tensor<4>& J, const Matrix& S, Matrix& D, bool print) const {
    // create all the necessary matrices, calculate ERI and initialize DIIS
    libint2::DIIS<Matrix> diis(this->diis.start, this->diis.keep);
    Matrix eps(H.rows(), 1), C(H.rows(), H.cols()); double E = 0;
    Tensor<4> ERI = J - 0.5 * J.shuffle(Array<4>{0, 3, 2, 1});

    // print the iteration header
    if (print) std::printf("ITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= maxiter; i++) {
        // define the contraction axes and start the timer
        Eigen::IndexPair<int> first(2, 0), second(3, 1);
        Timer::Timepoint start = Timer::Now();

        // calculate the fock matrix and extrapolate it
        Matrix F = H + toMatrix(ERI.contract(toTensor(D), Axes<2>{first, second}));
        Matrix e = S * D * F - F * D * S; if (i > 1) diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, S);
        C = solver.eigenvectors(), eps = solver.eigenvalues();
        Matrix Dp = D; double Ep = E;

        // calculate the new density and energy
        D = 2 * C.leftCols(system.nocc) * C.leftCols(system.nocc).transpose();
        E = 0.5 * D.cwiseProduct(H + F).sum();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, E, std::abs(E - Ep), (D - Dp).norm(), Timer::Format(Timer::Elapsed(start)).c_str(), i > this->diis.start + 1 ? "DIIS" : "");

        // finish if covergence reached
        if (std::abs(E - Ep) < thresh && (D - Dp).norm() < thresh) break;
    }

    // return the results
    return {C, eps, E};
}
