#include "../include/hartreefock.h"
#include "../include/molecule.h"
#include <libint2/diis.h>

typedef HartreeFockOptions Options;
typedef HartreeFockResult Result;

HartreeFock::HartreeFock(Options opt) : opt(opt) {
    // print method specification
    Printer::printMethod(opt);
}

HartreeFock::Result HartreeFock::scf(const Molecule& molecule) const {
    // initialize the result container
    Result result; auto times = &result.times;
    
    // calculate the necessary integrals
    auto start = Timer::now();
    Eigen::MatrixXd T = molecule.integral(libint2::Operator::kinetic);
    times->ints["T"]= Timer::elapsed(start), start = Timer::now();
    Eigen::MatrixXd V = molecule.integral(libint2::Operator::nuclear);
    times->ints["V"]= Timer::elapsed(start), start = Timer::now();
    Eigen::MatrixXd S = molecule.integral(libint2::Operator::overlap);
    times->ints["S"]= Timer::elapsed(start), start = Timer::now();
    double Vnn = molecule.getRepulsion(); Eigen::MatrixXd H = T + V;

    // save the integrals to the result container
    result.T = T, result.V = V, result.S = S, result.Vnn = Vnn;
    result.nocc = molecule.getElectronCount() / 2;

    // quess the initial density
    start = Timer::now();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H, S);
    Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
    Eigen::MatrixXd D =  C * C.transpose();
    times->guess.at(0) = Timer::elapsed(start);

    // compute the guess energy
    start = Timer::now();
    double E = D.cwiseProduct(2 * H).sum() + Vnn; start = Timer::now();
    times->guess.at(1) = Timer::elapsed(start);

    // print initial timings
    Printer::printInitialTimings(result);

    // initialize the DIIS algorithm and start the timer for the sCF cycle
    libint2::DIIS<Eigen::MatrixXd> diis(opt.diis.start - 1, opt.diis.keep, opt.diis.damp);
    start = Timer::now();

    // start the SCF cycle
    for (result.i = 1; result.i <= opt.maxiter; result.i++) {
        // compute the Fock matrix
        result.F = H + molecule.integral(libint2::Operator::coulomb, D);

        // compute error and extrapolate the Fock matrix
        Eigen::MatrixXd e = S * D * result.F - result.F * D * S, Fold = result.F;
        if (opt.diis.enabled) diis.extrapolate(result.F, e);

        // solve the Roothan equations and compute the density matrix from the result
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(result.F, S);
        Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
        result.D = C * C.transpose();

        // calculate the energy
        result.E = result.D.cwiseProduct(H + result.F).sum() + Vnn;

        // write the secondary results to the result container
        result.dD = std::abs(result.D.norm() - D.norm()), result.dE = std::abs(result.E - E);
        times->iters.push_back(Timer::elapsed(start)), start = Timer::now();

        // print the iteration
        Printer::printIteration(result, opt); D = result.D, E = result.E;

        // check for convergence
        if (result.dE < opt.thresh && result.dD < opt.thresh) { result.Eo = solver.eigenvalues(); break; }
        else if (result.i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    }
    // return the results
    return result;
}
