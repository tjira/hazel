#include "utility.h"

std::vector<std::vector<int>> Utility::Combinations(int n, int k) {
    // get the number of combinations and create resulting array
    std::vector<std::vector<int>> combs(Factorial(n) / (Factorial(k) * Factorial(n - k)));

    // create the bitmask that will get permuted
    std::string bitmask(k, 1); bitmask.resize(n, 0);
 
    // generate the combinations
    for (size_t i = 0; i < combs.size(); i++) {
        for (int j = 0; j < n; j++) {
            if (bitmask[j]) combs.at(i).push_back(j);
        }
        std::prev_permutation(bitmask.begin(), bitmask.end());
    }

    // return the result
    return combs;
}

std::vector<int> Utility::Common(const std::vector<int>& first, const std::vector<int>& second) {
    // create the container
    std::vector<int> common;

    // fill the container
    for (size_t i = 0; i < first.size(); i++) if (first.at(i) == second.at(i)) common.push_back(first.at(i));

    // return the container
    return common;
}

Matrix Utility::Repeat(const Matrix& M, int count, int axis) {
    Matrix N(axis == 0 ? count * M.rows() : M.rows(), axis == 1 ? count * M.cols() : M.cols());
    if (axis == 0) {
        for (int i = 0; i < M.rows(); i++) {
            for (int j = 0; j < count; j++) {
                N.row(i * count + j) = M.row(i);
            }
        }
    } else if (axis == 1) {
        for (int i = 0; i < M.cols(); i++) {
            for (int j = 0; j < count; j++) {
                N.col(i * count + j) = M.col(i);
            }
        }
    } else {
        throw std::runtime_error("UNKNOWN AXIS IN MATRIX REPEAT");
    }
    return N;
}

void Utility::SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn, const Vector& energy) {
    int size = 0;
    for (size_t i = 0; i < wfn.size(); i++) {
        size = std::max(size, (int)wfn.at(i).size());
    }
    std::ofstream file(fname);
    for (int i = 0; i < size; i++) {
        file << std::fixed << std::setprecision(14) << "#        r1         ";
        for (size_t j = 0; j < wfn.size(); j++) {
            file << "     state" << (j < 10 ? "0" : "") << j << ".real         " << "state" << (j < 10 ? "0" : "") << j << ".imag    ";
        }
        if (!i) {
            file << " E=[";
            for (int j = 0; j < energy.size(); j++) file << energy(j) << (j < energy.size() - 1 ? ", " : "]");
        }
        for (long int j = 0; j < r.size(); j++) {
            file << (!j ? "\n" : "") << std::setw(20) << r(j).real();
            for (size_t k = 0; k < wfn.size(); k++) {
                file << " " << std::setw(20) << wfn.at(k).at(std::min(i, (int)wfn.at(k).size() - 1))(j).real()
                << " " << std::setw(20) << wfn.at(k).at(std::min(i, (int)wfn.at(k).size() - 1))(j).imag();
            }
            file << "\n";
        }
    }
}

std::vector<int> Utility::Unique(const std::vector<int>& first, const std::vector<int>& second) {
    // create the container
    std::vector<int> unique;

    // fill the container
    for (size_t i = 0; i < second.size(); i++) if (!VectorContains(first, second.at(i))) unique.push_back(second.at(i));
    for (size_t i = 0; i < first.size(); i++) if (!VectorContains(second, first.at(i))) unique.push_back(first.at(i));

    // return the container
    return unique;
}
