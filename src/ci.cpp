#include "ci.h"

Data CI::cid(bool) const {
    // calculate ERI and define output
    Tensor<4> KMO = data.intsmo.J.shuffle(Array<4>{0, 3, 2, 1});
    Data output = data; int nocc = data.system.electrons / 2;

    // generate extications
    std::vector<std::pair<int, int>> doubles;
    for (int i = 0; i < nocc; i++) {
        for (int j = nocc; j < data.system.shells.nbf(); j++) {
            doubles.push_back({i, j});
        }
    }

    // create the CI Hamiltonian an define the contraction axes
    output.ci.H = Matrix::Zero(doubles.size() + 1, doubles.size() + 1);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // fill the first CI matrix element
    output.ci.H(0, 0) = data.hf.E - Integral::Repulsion(data.system);

    // fill the double in CI Hamiltonian
    for (size_t i = 0; i < doubles.size(); i++) {
        // extract extication level and element index
        int a = doubles.at(i).first, b = doubles.at(i).second;

        // assign the diagonal element
        output.ci.H(i + 1, i + 1) = output.ci.H(0, 0) + 2 * (data.hf.eps(b) - data.hf.eps(a)) + data.intsmo.J(a, a, a, a) + data.intsmo.J(b, b, b, b) - 4 * data.intsmo.J(b, b, a, a) + 2 * KMO(b, b, a, a);

        // fill the matrix elements with HF determinant
        output.ci.H(0, i + 1) = KMO(doubles.at(i).first, doubles.at(i).first, doubles.at(i).second, doubles.at(i).second);
        output.ci.H(i + 1, 0) = KMO(doubles.at(i).second, doubles.at(i).second, doubles.at(i).first, doubles.at(i).first);

        // fill double-double matrix elements
        for (size_t j = 0; j < i; j++) {
            output.ci.H(j + 1, i + 1) = KMO(doubles.at(j).second, doubles.at(j).second, doubles.at(i).second, doubles.at(i).second);
            output.ci.H(i + 1, j + 1) = KMO(doubles.at(j).second, doubles.at(j).second, doubles.at(i).second, doubles.at(i).second);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    Eigen::SelfAdjointEigenSolver<Matrix> solver(output.ci.H);

    // extract the excitation energies and expansion coefs
    output.ci.C = solver.eigenvectors(), output.ci.eig = solver.eigenvalues().array() + Integral::Repulsion(data.system);

    // extract correlation energy
    output.ci.Ecorr = output.ci.eig(0) - data.hf.E;

    // return the result
    return output;
}

Data CI::cis(bool) const {
    // define the output and occupation number
    Data output = data; int nocc = data.system.electrons / 2;

    // generate extications
    std::vector<std::pair<int, int>> singles;
    for (int i = 0; i < nocc; i++) {
        for (int j = nocc; j < data.system.shells.nbf(); j++) {
            singles.push_back({i, j});
        }
    }

    // reverse vector
    std:reverse(singles.begin(), singles.end());

    // create the CI Hamiltonian an define the contraction axes
    output.ci.H = Matrix::Zero(singles.size() + 1, singles.size() + 1);
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // fill the first CI matrix element
    output.ci.H(0, 0) = data.hf.E - Integral::Repulsion(data.system);

    // fill the singlet singles in CI Hamiltonian
    for (size_t i = 0; i < singles.size(); i++) {
        for (size_t j = 0; j < singles.size(); j++) {
            // extract extication levels
            int a = singles.at(i).first, r = singles.at(i).second, b = singles.at(j).first, s = singles.at(j).second;

            // assign the elements
            if (i == j) output.ci.H(i + 1, j + 1) = output.ci.H(0, 0) + (data.hf.eps(r) - data.hf.eps(a)) - data.intsmo.J(r, s, b, a) + data.intsmo.J(r, a, b, s);
            else output.ci.H(i + 1, j + 1) = data.intsmo.J(r, a, b, s) - data.intsmo.J(r, s, b, a);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian
    Eigen::SelfAdjointEigenSolver<Matrix> solver(output.ci.H);

    // extract the excitation energies and expansion coefs
    output.ci.C = solver.eigenvectors(), output.ci.eig = solver.eigenvalues().array() + Integral::Repulsion(data.system);

    // extract correlation energy
    output.ci.Ecorr = output.ci.eig(0) - data.hf.E;

    // return the result
    return output;
}
