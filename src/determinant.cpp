#include "determinant.h"

Determinant::Determinant(int norb, int nocca, int noccb) : a(nocca), b(noccb), norb(norb) {
    std::iota(a.begin(), a.end(), 0), std::iota(b.begin(), b.end(), 0);
}

Determinant::Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b) : a(a), b(b), norb(norb) {}

std::ostream& operator<<(std::ostream& os, const Determinant& det) {
    os << "|(";
    for (size_t i = 0; i < det.a.size(); i++) {
        if (i == det.a.size() - 1) os << det.a.at(i);
        else os << det.a.at(i) << ",";
    }
    os << "),(";
    for (size_t i = 0; i < det.b.size(); i++) {
        if (i == det.b.size() - 1) os << det.b.at(i);
        else os << det.b.at(i) << ",";
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

double Determinant::hamilton(const Determinant& det, const Matrix& Hms, const Tensor<4>& Jms) const {
    // create all needed variables
    std::vector<std::vector<int>> exc; int swaps = 0, diff = 0; double elem = 0;
    std::vector<int> SOorig, SOdet, common; Determinant orig(*this);
    
    // align the alpha electron in determinants
    for (size_t i = 0; i < orig.a.size(); i++) {
        if (orig.a.at(i) != det.a.at(i)) {
            for (size_t j = i + 1; j < a.size(); j++) {
                if (orig.a.at(i) == det.a.at(j)) {
                    std::swap(orig.a.at(i), orig.a.at(j)), swaps++, i = -1; break;
                }
            }
        }
    }

    // align the beta electron in determinants
    for (size_t i = 0; i < orig.b.size(); i++) {
        if (orig.b.at(i) != det.b.at(i)) {
            for (size_t j = i + 1; j < b.size(); j++) {
                if (orig.b.at(i) == det.b.at(j)) {
                    std::swap(orig.b.at(i), orig.b.at(j)), swaps++, i = -1; break;;
                }
            }
        }
    }

    // fill the spinorbitals of the original determinant
    for (int i = 0; i < orig.b.size(); i++) SOorig.push_back(2 * orig.b.at(i) + 1);
    for (int i = 0; i < orig.a.size(); i++) SOorig.push_back(2 * orig.a.at(i));

    // fill the spinorbitals of the second determinant
    for (int i = 0; i < det.b.size(); i++) SOdet.push_back(2 * det.b.at(i) + 1);
    for (int i = 0; i < det.a.size(); i++) SOdet.push_back(2 * det.a.at(i));

    // fill in the common orbitals
    for (int i = 0; i < SOorig.size(); i++) if (SOorig.at(i) == SOdet.at(i)) common.push_back(SOorig.at(i));

    // get the number of different spinorbitals
    for (int i = 0; i < orig.a.size(); i++) if (orig.a.at(i) != det.a.at(i)) diff++;
    for (int i = 0; i < orig.b.size(); i++) if (orig.b.at(i) != det.b.at(i)) diff++;

    // extract excitations and make sure to convert spatial orbital to spinorbital
    for (int i = 0; i < orig.b.size(); i++) if (orig.b.at(i) > det.b.at(i)) exc.push_back({2 * det.b.at(i) + 1, 2 * orig.b.at(i) + 1});
    for (int i = 0; i < orig.b.size(); i++) if (orig.b.at(i) < det.b.at(i)) exc.push_back({2 * orig.b.at(i) + 1, 2 * det.b.at(i) + 1});
    for (int i = 0; i < orig.a.size(); i++) if (orig.a.at(i) > det.a.at(i)) exc.push_back({2 * det.a.at(i), 2 * orig.a.at(i)});
    for (int i = 0; i < orig.a.size(); i++) if (orig.a.at(i) < det.a.at(i)) exc.push_back({2 * orig.a.at(i), 2 * det.a.at(i)});

    // assign the matrix element
    if (diff == 0) {
        for (int so : SOorig) elem += Hms(so, so);
        for (size_t i = 0; i < SOorig.size(); i++) {
            for (size_t j = i + 1; j < SOorig.size(); j++) {
                elem += Jms(SOorig.at(i), SOorig.at(j), SOorig.at(i), SOorig.at(j));
            }
        }
    } else if (diff == 1) {
        elem += Hms(exc.at(0).at(0), exc.at(0).at(1));
        for (int so : common) {
            elem += Jms(exc.at(0).at(0), so, exc.at(0).at(1), so);
        }
    } else if (diff == 2) {
        elem = Jms(exc.at(0).at(1), exc.at(1).at(1), exc.at(0).at(0), exc.at(1).at(0));
    }

    // return the matrix element
    return std::pow(-1, swaps) * elem;
}
