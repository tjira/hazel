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
    Eigen::MatrixXd T = molecule.kinetic(), V = molecule.nuclear();
    Eigen::MatrixXd S = molecule.overlap(), H = T + V;
    double Vnn = molecule.nuclearRepulsion();
    times->ints = Timer::elapsed(start);

    result.T = T, result.V = V, result.S = S, result.Vnn = Vnn;
    int nocc = molecule.nel() / 2;

    start = Timer::now();
    Eigen::MatrixXd D = computeDensity(S, H, molecule.nel() / 2);
    result.Es.push_back(computeEnergy(H, H, D) + Vnn);
    result.Fs.push_back(H), result.Ds.push_back(D);
    result.DNs.push_back(D.norm());
    times->guess = Timer::elapsed(start);

    start = Timer::now();
    for (int i = 1; i <= opt.maxiter; i++) {
        Eigen::MatrixXd F = H + molecule.coulomb(D);
        D = computeDensity(S, F, nocc);

        result.Es.push_back(computeEnergy(H, F, D) + Vnn);
        result.Fs.push_back(F), result.Ds.push_back(D);
        result.DNs.push_back(D.norm());

        if (opt.damp > 0) D = (1 - opt.damp) * D + opt.damp * result.Ds.at(i - 1);
        times->iters.push_back(Timer::elapsed(start)), start = Timer::now();

        bool converged = checkConvergence(result, i);
        Printer::printIteration(result, i, converged);
        if (converged) break;
    }
    Printer::printResult(result); return result;
}

bool HartreeFock::checkConvergence(const Result& result, int i) const {
    if (std::abs(result.Es.at(i - 1) - result.Es.at(i)) < opt.thresh && std::abs(result.DNs.at(i - 1) - result.DNs.at(i)) < opt.thresh) return true;
    else if (i == opt.maxiter) std::cerr << "Algorithm did not converge." << std::endl;
    return false;
}

Eigen::MatrixXd HartreeFock::computeDensity(Eigen::MatrixXd S, Eigen::MatrixXd F, int nocc) const {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(F, S);
    Eigen::MatrixXd C = solver.eigenvectors().leftCols(nocc);
    return C * C.transpose();
}

double HartreeFock::computeEnergy(Eigen::ArrayXXd H, Eigen::ArrayXXd F, Eigen::ArrayXXd D) const {
    return (D * (H + F)).sum();
};

void HartreeFock::logIteration(int i, const Result& result) const {
    std::cout << std::setw(4) << i + 1 << " " << result.Es.at(i) - result.Vnn << " " << result.Es.at(i);
    if (i > 0) std::cout << " " << std::abs(result.Es.at(i - 1) - result.Es.at(i));
    std::cout << std::endl;
}
