#include "ci.h"

Data CI::cid(bool) const {
    // calculate ERI and K in MO basis
    Tensor<4> ERI = data.ints.J - 0.5 * data.ints.J.shuffle(Array<4>{0, 3, 2, 1});
    Tensor<4> KMO = data.intsmo.J.shuffle(Array<4>{0, 3, 2, 1});
    Data output = data; int nocc = data.system.electrons / 2;

    // generate extications
    std::vector<std::pair<int, int>> excs;
    for (int i = 0; i < nocc; i++) {
        for (int j = nocc; j < data.system.shells.nbf(); j++) {
            excs.push_back({i, j});
        }
    }

    // create the CI Hamiltonian an define the contraction axes
    output.ci.H = Matrix::Zero(excs.size() + 1, excs.size() + 1);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // fill the first CI matrix element
    output.ci.H(0, 0) = output.hf.E - Integral::Repulsion(data.system);

    // fill the CI Hamiltonian
    for (size_t i = 0; i < excs.size(); i++) {
        // create excited coefficient matrix and density
        Matrix CEXC = data.hf.C; CEXC.col(excs.at(i).first).swap(CEXC.col(excs.at(i).second));
        Matrix DEXC = 2 * CEXC.leftCols(nocc) * CEXC.leftCols(nocc).transpose();

        // calculate the fock matrix end excitated energy
        Matrix FEXC = data.ints.T + data.ints.V + toMatrix(ERI.contract(toTensor(DEXC), Axes<2>{first, second}));
        output.ci.H(i + 1, i + 1) = 0.5 * DEXC.cwiseProduct(data.ints.T + data.ints.V + FEXC).sum();

        // fill the matrix elements with HF determinant
        output.ci.H(0, i + 1) = KMO(excs.at(i).first, excs.at(i).first, excs.at(i).second, excs.at(i).second);
        output.ci.H(i + 1, 0) = KMO(excs.at(i).second, excs.at(i).second, excs.at(i).first, excs.at(i).first);

        // fill double-double matrix elements
        for (size_t j = 0; j < i; j++) {
            output.ci.H(j + 1, i + 1) = KMO(excs.at(j).second, excs.at(j).second, excs.at(i).second, excs.at(i).second);
            output.ci.H(i + 1, j + 1) = output.ci.H(j + 1, i + 1);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    Eigen::SelfAdjointEigenSolver<Matrix> solver(output.ci.H);

    // extract the excitation energies and expansion coefs
    output.ci.C = solver.eigenvectors(), output.ci.eig = solver.eigenvalues();

    // extract correlation energy
    output.ci.Ecorr = output.ci.eig(0) - data.hf.E + Integral::Repulsion(data.system);

    // return the result
    return output;
}
