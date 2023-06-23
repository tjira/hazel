#include "ci.h"

CI::CI(const System& system, double Eel) : Eel(Eel), nocc(system.electrons / 2), nbf(system.shells.nbf()) {}

std::tuple<Matrix, Matrix, Vector> CI::cid(const Integrals& ints, const Tensor<4>& JMO, const Matrix& C) const {
    // calculate ERI and K in MO basis
    Tensor<4> ERI = ints.J - 0.5 * ints.J.shuffle(Array<4>{0, 3, 2, 1});
    Tensor<4> KMO = JMO.shuffle(Array<4>{0, 3, 2, 1});

    // generate extications
    std::vector<std::pair<int, int>> excs;
    for (int i = 0; i < nocc; i++) {
        for (int j = nocc; j < nbf; j++) {
            excs.push_back({i, j});
        }
    }

    // create the CI Hamiltonian an define the contraction axes
    Matrix CIH(excs.size() + 1, excs.size() + 1); CIH(0, 0) = Eel;
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // fill the CI Hamiltonian
    for (size_t i = 0; i < excs.size(); i++) {
        // create excited coefficient matrix and density
        Matrix CEXC = C; CEXC.col(excs.at(i).first).swap(CEXC.col(excs.at(i).second));
        Matrix DEXC = 2 * CEXC.leftCols(nocc) * CEXC.leftCols(nocc).transpose();

        // calculate the fock matrix end excitated energy
        Matrix FEXC = ints.T + ints.V + toMatrix(ERI.contract(toTensor(DEXC), Axes<2>{first, second}));
        CIH(i + 1, i + 1) = 0.5 * DEXC.cwiseProduct(ints.T + ints.V + FEXC).sum();

        // fill the matrix elements with HF determinant
        CIH(0, i + 1) = KMO(excs.at(i).first, excs.at(i).first, excs.at(i).second, excs.at(i).second);
        CIH(i + 1, 0) = KMO(excs.at(i).second, excs.at(i).second, excs.at(i).first, excs.at(i).first);

        // fill double-double matrix elements
        for (size_t j = 0; j < i; j++) {
            CIH(j + 1, i + 1) = KMO(excs.at(j).second, excs.at(j).second, excs.at(i).second, excs.at(i).second);
            CIH(i + 1, j + 1) = CIH(j + 1, i + 1);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    Eigen::SelfAdjointEigenSolver<Matrix> solver(CIH);

    // return the energy
    return {CIH, solver.eigenvectors(), solver.eigenvalues()};
}
