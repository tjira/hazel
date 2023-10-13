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

std::tuple<Determinant, int> Determinant::align(const Determinant& second) const {
    // define the temporary determinant and swaps variable
    Determinant first(*this); int swaps = 0;
    
    // align the alpha electrons
    for (size_t i = 0; i < first.a.size(); i++) {
        if (first.a.at(i) != second.a.at(i)) {
            for (size_t j = 0; j < a.size(); j++) {
                if (first.a.at(i) == second.a.at(j)) {
                    std::swap(first.a.at(i), first.a.at(j)), swaps++;
                }
            }
        }
    }

    // align the beta electrons
    for (size_t i = 0; i < first.b.size(); i++) {
        if (first.b.at(i) != second.b.at(i)) {
            for (size_t j = 0; j < b.size(); j++) {
                if (first.b.at(i) == second.b.at(j)) {
                    std::swap(first.b.at(i), first.b.at(j)), swaps++;
                }
            }
        }
    }

    // return the aligned determinant
    return std::tuple<Determinant, int>{first, swaps};
}

int Determinant::differences(const Determinant& second) const {
    // define the count
    int count = 0;

    // count the differences
    for (size_t i = 0; i < a.size(); i++) if (a.at(i) != second.a.at(i)) count++;
    for (size_t i = 0; i < b.size(); i++) if (b.at(i) != second.b.at(i)) count++;

    // return the count
    return count;
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
    // align this determinant, number of swaps and differences and the matrix element
    auto[first, swaps] = align(second); int diff = first.differences(second); double elem = 0;
    std::vector<int> firstso = first.spinorbitals(), secondso = second.spinorbitals();

    // get the common and unique  spinorbitals
    std::vector<int> common = Utility::Common(firstso, secondso), unique = Utility::Unique(firstso, secondso);

    // assign the matrix element
    if (diff == 0) {
        for (int so : firstso) elem += Hms(so, so);
        for (size_t i = 0; i < firstso.size(); i++) {
            for (size_t j = i + 1; j < firstso.size(); j++) {
                elem += Jms(firstso.at(i), firstso.at(i), firstso.at(j), firstso.at(j)) - Jms(firstso.at(i), firstso.at(j), firstso.at(j), firstso.at(i));
            }
        }
    } else if (diff == 1) {
        elem += Hms(unique.at(0), unique.at(1));
        for (int so : common) {
            elem += Jms(unique.at(0), unique.at(1), so, so) - Jms(unique.at(0), so, unique.at(1), so);
        }
    } else if (diff == 2) {
        elem = Jms(unique.at(0), unique.at(2), unique.at(1), unique.at(3)) - Jms(unique.at(0), unique.at(3), unique.at(2), unique.at(1));
    }

    // return diff;
    return std::pow(-1, swaps) * elem;
}

std::vector<int> Determinant::spinorbitals() const {
    // create the vector to hold all spinorbitals
    std::vector<int> spinorbitals(a.size() + b.size());

    // fill the spinorbitals vector
    for (size_t i = 0; i < a.size(); i++) spinorbitals.at(b.size() + i) = 2 * a.at(i);
    for (size_t i = 0; i < b.size(); i++) spinorbitals.at(i) = 2 * b.at(i) + 1;

    // return the spinorbital vector
    return spinorbitals;
}
