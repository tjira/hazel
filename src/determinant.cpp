#include "determinant.h"

Determinant::Determinant(int norb, int nocca, int noccb) : a(nocca), b(noccb), norb(norb) {
    std::iota(a.begin(), a.end(), 0), std::iota(b.begin(), b.end(), 0);
}

Determinant::Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b) : a(a), b(b), norb(norb) {}

std::ostream& operator<<(std::ostream& os, const Determinant& second) {
    os << "|(";
    for (size_t i = 0; i < second.a.size(); i++) {
        if (i == second.a.size() - 1) os << second.a.at(i);
        else os << second.a.at(i) << ",";
    }
    os << "),(";
    for (size_t i = 0; i < second.b.size(); i++) {
        if (i == second.b.size() - 1) os << second.b.at(i);
        else os << second.b.at(i) << ",";
    }
    os << ")>";
    return os;
}

std::vector<Determinant> Determinant::full() const {
    std::vector<Determinant> full;
    for (const std::vector<int>& alpha : Utility::Combinations(norb, a.size())) {
        for (const std::vector<int>& beta : Utility::Combinations(norb, b.size())) {
            full.push_back(Determinant(norb, alpha, beta));
        }
    }
    return full;
}

double Determinant::hamilton(const Determinant& second, const Matrix& Hms, const Tensor<4>& Jms) const {
    // create all needed variables
    int swaps = 0, diff = 0; double elem = 0; Determinant first(*this);
    std::vector<int> firstso, secondso, common, unique;
    
    // align the alpha electron in seterminants
    for (size_t i = 0; i < first.a.size(); i++) {
        if (first.a.at(i) != second.a.at(i)) {
            for (size_t j = 0; j < a.size(); j++) {
                if (first.a.at(i) == second.a.at(j)) {
                    std::swap(first.a.at(i), first.a.at(j)), swaps++;
                }
            }
        }
    }

    // align the beta electron in determinants
    for (size_t i = 0; i < first.b.size(); i++) {
        if (first.b.at(i) != second.b.at(i)) {
            for (size_t j = 0; j < b.size(); j++) {
                if (first.b.at(i) == second.b.at(j)) {
                    std::swap(first.b.at(i), first.b.at(j)), swaps++;
                }
            }
        }
    }

    // fill the spinorbitals of the firstinal determinant
    for (size_t i = 0; i < first.b.size(); i++) firstso.push_back(2 * first.b.at(i) + 1);
    for (size_t i = 0; i < first.a.size(); i++) firstso.push_back(2 * first.a.at(i));

    // fill the spinorbitals of the second determinant
    for (size_t i = 0; i < second.b.size(); i++) secondso.push_back(2 * second.b.at(i) + 1);
    for (size_t i = 0; i < second.a.size(); i++) secondso.push_back(2 * second.a.at(i));

    // fill in the common orbitals
    for (size_t i = 0; i < firstso.size(); i++) if (firstso.at(i) == secondso.at(i)) common.push_back(firstso.at(i));

    // get the number of different spinorbitals
    for (size_t i = 0; i < first.a.size(); i++) if (first.a.at(i) != second.a.at(i)) diff++;
    for (size_t i = 0; i < first.b.size(); i++) if (first.b.at(i) != second.b.at(i)) diff++;

    // et number of unique orbitals
    for (size_t i = 0; i < firstso.size(); i++) {
        if (std::find(secondso.begin(), secondso.end(), firstso.at(i)) == secondso.end()) {
            unique.push_back(firstso.at(i));
        }
    }
    for (size_t i = 0; i < secondso.size(); i++) {
        if (std::find(firstso.begin(), firstso.end(), secondso.at(i)) == firstso.end()) {
            unique.push_back(secondso.at(i));
        }
    }

    // assign the matrix element
    if (diff == 0) {
        for (int so : firstso) elem += Hms(so, so);
        for (size_t i = 0; i < firstso.size(); i++) {
            for (size_t j = i + 1; j < firstso.size(); j++) {
                elem += Jms(firstso.at(i), firstso.at(j), firstso.at(i), firstso.at(j));
            }
        }
    } else if (diff == 1) {
        elem += Hms(unique.at(0), unique.at(1));
        for (int so : common) {
            elem += Jms(unique.at(0), so, unique.at(1), so);
        }
    } else if (diff == 2) {
        elem = Jms(unique.at(0), unique.at(1), unique.at(2), unique.at(3));
    }

    // return diff;
    return std::pow(-1, swaps) * elem;
}
