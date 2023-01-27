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

    result.Es.push_back((2 * D.array() * H.array()).sum() + Vnn);
    result.Fs.push_back(H), result.Ds.push_back(D);

    start = Timer::now();
    for (result.iters = 1; result.iters <= opt.maxiter; result.iters++) {
        Eigen::MatrixXd F = H + molecule.integral<2>(libint2::Operator::coulomb, D);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(F, S);
        Eigen::MatrixXd C = solver.eigenvectors().leftCols(result.nocc);
        D = C * C.transpose();

        result.Es.push_back((D.array() * (H + F).array()).sum() + Vnn);
        result.Fs.push_back(F), result.Ds.push_back(D);
        result.Eo = solver.eigenvalues();

        if (opt.damp > 0) D = (1 - opt.damp) * D + opt.damp * result.Ds.at(result.iters - 1);
        times->iters.push_back(Timer::elapsed(start)), start = Timer::now();

        double dD = std::abs(result.Ds.at(result.iters).norm() - result.Ds.at(result.iters - 1).norm());
        double dE = std::abs(result.Es.at(result.iters) - result.Es.at(result.iters - 1));

        bool converged = dE < opt.thresh && dD < opt.thresh, lastiter = result.iters == opt.maxiter;
        Printer::printIteration(result, converged || lastiter); if (converged) break;
        if (lastiter) std::cerr << "Algorithm did not converge." << std::endl;
    }
    return result;
}
