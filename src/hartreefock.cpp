#include "../include/hartreefock.h"
#include "../include/molecule.h"
#include <boost/format.hpp>
#include <libint2/diis.h>

typedef HartreeFockOptions Options;
typedef HartreeFockResult Result;

HartreeFock::HartreeFock(Options opt) : opt(opt) {
    // print method specification
    std::cout << "HARTREE-FOCK" << std::endl;
    std::cout << boost::format("MAXITER: %i, THRESH: %.2e, DIIS: [ENABLED: %i, START: %i, KEEP: %i, DAMP: %.2f]") % opt.maxiter
    % opt.thresh % opt.diis.enabled % opt.diis.start % opt.diis.keep % opt.diis.damp << std::endl << std::endl;
}

HartreeFock::Result HartreeFock::scf(const Molecule& molecule) const {
    // print initial info
    std::cout << "SCF CYCLE" << std::endl;
    std::cout << "ITER        E [Eh]           dE       dD        TIME" << std::endl;

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
    double Vnn = molecule.getNuclearRepulsion(); Eigen::MatrixXd H = T + V;

    // save the integrals to the result container
    result.T = T, result.V = V, result.S = S, result.Vnn = Vnn;
    result.nocc = molecule.getElectronCount() / 2;

    // quess the initial density
    start = Timer::now();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H, S);
    Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
    Eigen::MatrixXd D =  2 * C * C.transpose();
    times->guess.at(0) = Timer::elapsed(start);

    // compute the guess energy
    start = Timer::now();
    double E = D.cwiseProduct(2 * H).sum() + Vnn; start = Timer::now();
    times->guess.at(1) = Timer::elapsed(start);

    // initialize the DIIS algorithm and start the timer for the SCF cycle
    libint2::DIIS<Eigen::MatrixXd> diis(opt.diis.start - 1, opt.diis.keep, opt.diis.damp);
    start = Timer::now();

    // start the SCF cycle
    for (result.i = 1; result.i <= opt.maxiter; result.i++) {
        // compute the Fock matrix
        result.F = H + 0.5 * molecule.integral(libint2::Operator::coulomb, D);

        // compute error and extrapolate the Fock matrix
        Eigen::MatrixXd e = S * D * result.F - result.F * D * S, Fold = result.F;
        if (opt.diis.enabled) diis.extrapolate(result.F, e);

        // solve the Roothan equations and compute the density matrix from the result
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(result.F, S);
        Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
        result.C = C, result.D = 2 * C * C.transpose();

        // calculate the energy
        result.E = 0.5 * result.D.cwiseProduct(H + result.F).sum() + Vnn;

        // write the secondary results to the result container
        result.dD = std::abs(result.D.norm() - D.norm()), result.dE = std::abs(result.E - E);
        times->iters.push_back(Timer::elapsed(start)), start = Timer::now();

        // print the iteration
        std::cout << boost::format("%4i %20.14f %.2e %.2e %s") % result.i % result.E % result.dE % result.dD
        % Timer::format(times->iters.at(result.i - 1)) << std::endl;

        // overwrite matrices
        D = result.D, E = result.E;

        // check for convergence
        if (result.dE < opt.thresh && result.dD < opt.thresh) { result.eps = solver.eigenvalues(); break; }
        else if (result.i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    }
    // print blank line and return the results
    std::cout << std::endl; return result;
}
