#include "../include/hartreefock.h"
#include "../include/molecule.h"

typedef HartreeFockOptions Options;
typedef HartreeFockResult Result;

HartreeFock::HartreeFock(Options opt) : opt(opt) {
    Printer::printMethod(opt);
}

HartreeFock::Result HartreeFock::scf(const Molecule& molecule) const {
    Result result; auto times = &result.times;
    
    auto start = Timer::now();
    Eigen::MatrixXd T = molecule.integral<1>(libint2::Operator::kinetic);
    Eigen::MatrixXd V = molecule.integral<1>(libint2::Operator::nuclear);
    Eigen::MatrixXd S = molecule.integral<1>(libint2::Operator::overlap);
    double Vnn = molecule.nuclearRepulsion(); Eigen::MatrixXd H = T + V;
    times->ints = Timer::elapsed(start);

    result.T = T, result.V = V, result.S = S, result.Vnn = Vnn;
    result.nocc = molecule.nel() / 2;

    start = Timer::now();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H, S);
    Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
    Eigen::MatrixXd D =  C * C.transpose();
    times->guess = Timer::elapsed(start);

    libint2::DIIS<Eigen::MatrixXd> diis(opt.diis.start - 1, opt.diis.keep, opt.diis.damp);
    double E = D.cwiseProduct(2 * H).sum() + Vnn; start = Timer::now();

    for (result.i = 1; result.i <= opt.maxiter; result.i++) {
        result.F = H + molecule.integral<2>(libint2::Operator::coulomb, D);
        Eigen::MatrixXd e = S * D * result.F - result.F * D * S, Fold = result.F;
        if (opt.diis.start > 0) diis.extrapolate(result.F, e);

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(result.F, S);
        Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
        result.D = C * C.transpose();

        if (opt.damp > 0) result.D = (1 - opt.damp) * result.D + opt.damp * D;
        result.E = result.D.cwiseProduct(H + result.F).sum() + Vnn;

        result.dD = std::abs(result.D.norm() - D.norm()), result.dE = std::abs(result.E - E);
        times->iters.push_back(Timer::elapsed(start)), start = Timer::now();
        Printer::printIteration(result, opt); D = result.D, E = result.E;

        if (result.dE < opt.thresh && result.dD < opt.thresh) { result.Eo = solver.eigenvalues(); break; }
        else if (result.i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    }
    return result;
}
