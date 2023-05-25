#include "roothaan.h"

Roothaan::Roothaan(int maxiter, double thresh) : thresh(thresh), maxiter(maxiter) {}

std::tuple<Matrix, Matrix, double> Roothaan::scf(const Matrix& H, const Tensor<4>& J, const Matrix& S, Matrix D, int nocc, bool print) const {
    // print the method header
    if (print) std::cout << "\n" + std::string(WIDTH, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n";
    if (print) std::cout << std::format("MAXITER: {:d}, THRESH: {:.2e}, DIIS: [START: {:d}, KEEP: {:d}]\n", maxiter, thresh, 3, 5);
    if (print) std::cout << std::string(WIDTH, '-') + "\n\n";

    // create all the necessary matrices and vectors
    Tensor<4> ERI = J - 0.5 * J.t({0, 3, 2, 1});
    Matrix eps(H.rows(), 1); double E = 0;
    Matrix C(H.rows(), H.cols());

    // initialize DIIS
    libint2::DIIS<EigenMat> diis(3, 5);

    // print the iteration header
    if (print) std::cout << std::format("ITER {:^20s} {:^8s} {:^8s} {:^12s}", "Eel [Eh]", "|dE|", "|dD|", "TIME") << std::endl;

    // perform the scf loop
    for (int i = 1; i <= maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // calculate the fock matrix and the error vector
        Matrix F = H + ERI.contract(D, {Ind{2, 0}, Ind{3, 1}});
        Matrix e = S.dot(D).dot(F) - F.dot(D).dot(S);

        // extrapolate the fock matrix and solve the roothaan equations
        if (i > 1) diis.extrapolate(F.raw(), e.raw());
        std::tie(C, eps) = F.eigensolve(S);
        Matrix Dp = D; double Ep = E;

        // calculate the new density and energy
        D = 2 * C.left(nocc).dot(C.left(nocc).t());
        E = 0.5 * (D * (H + F)).sum();

        // print the iteration info line
        if (print) std::cout << std::format("{:4d} {:20.14f} {:.2e} {:.2e} {:12s}", i, E, std::abs(E - Ep), (D - Dp).norm(), Timer::Format(Timer::Elapsed(start))) << std::endl;

        // finish if covergence reached
        if (std::abs(E - Ep) < thresh && (D - Dp).norm() < thresh) break;
    }

    // return the results
    return {C, eps, E};
}
