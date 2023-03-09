#include "../include/hartreefock.h"
#include <Eigen/src/Core/Matrix.h>
#include <libint2/atom.h>
#include <valarray>

System::System(std::string filename, std::string basis) : basis(basis) {
    std::ifstream file(filename); atoms = libint2::read_dotxyz(file), shells = libint2::BasisSet(basis, atoms, true);
    /* std::copy(shells.begin(), shells.end(), std::ostream_iterator<libint2::Shell>(std::cout, "\n")); */
}


// algorithm from: https://www.cup.uni-muenchen.de/ch/compchem/pop/mull1.html
MullikenResult System::mulliken(Eigen::MatrixXd D) const {
    Eigen::MatrixXd S = integralSingle(libint2::Operator::overlap);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(atoms.size());
    Eigen::MatrixXd DS = D.cwiseProduct(S);
    for (size_t i = 0, j = 0; i < shells.size(); i++) {
        for (size_t k = 0; k < shells.at(i).size(); k++) {
            q(shells.shell2atom(atoms).at(i)) -= DS.colwise().sum()(j + k);
        }
        j += shells.at(i).size();
    }
    for (size_t i = 0; i < atoms.size(); i++) {
        q(i) += atoms.at(i).atomic_number;
    }
    return { DS, q };
}

/*
Computes Drs * (2 * (pq|rs) - (pr|sq))).
*/
Eigen::MatrixXd System::integralCoulomb(Eigen::MatrixXd D) const {
    libint2::Engine engine(libint2::Operator::coulomb, shells.max_nprim(), shells.max_l());
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(shells.nbf(), shells.nbf());
    for(size_t i = 0; i < shells.size(); i++) {
        const auto& result = engine.results(); auto sh2bf = shells.shell2bf();
        for(size_t j = 0; j <= i; j++) {
            for(size_t k = 0; k <= i; k++) {
                for(size_t l = 0; l <= (i == k ? j : k); l++) {
                    engine.compute(shells[i], shells[j], shells[k], shells[l]); if (result[0] == nullptr) continue;
                    double degeneracy = (i == j ? 1 : 2) * (k == l ? 1 : 2) * (i == k ? (j == l ? 1 : 2) : 2);
                    for(size_t m = 0, q = 0; m < shells[i].size(); m++) {
                        for(size_t n = 0; n < shells[j].size(); n++) {
                            for(size_t o = 0; o < shells[k].size(); o++) {
                                for(size_t p = 0; p < shells[l].size(); p++, q++) {
                                    size_t bf1 = m + sh2bf[i], bf2 = n + sh2bf[j];
                                    size_t bf3 = o + sh2bf[k], bf4 = p + sh2bf[l];
                                    matrix(bf1, bf3) -= 0.25 * D(bf2, bf4) * result[0][q] * degeneracy;
                                    matrix(bf2, bf4) -= 0.25 * D(bf1, bf3) * result[0][q] * degeneracy;
                                    matrix(bf1, bf4) -= 0.25 * D(bf2, bf3) * result[0][q] * degeneracy;
                                    matrix(bf2, bf3) -= 0.25 * D(bf1, bf4) * result[0][q] * degeneracy;
                                    matrix(bf1, bf2) += D(bf3, bf4) * result[0][q] * degeneracy;
                                    matrix(bf3, bf4) += D(bf1, bf2) * result[0][q] * degeneracy;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0.5 * (matrix + matrix.transpose());
};

Eigen::MatrixXd System::integralSingle(libint2::Operator op) const {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(shells.nbf(), shells.nbf());
    libint2::Engine engine(op, shells.max_nprim(), shells.max_l());
    if (op == libint2::Operator::nuclear) {
        engine.set_params(libint2::make_point_charges(atoms));
    }
    for (size_t i = 0; i < shells.size(); i++) {
        const auto& result = engine.results(); auto sh2bf = shells.shell2bf();
        for (size_t j = 0; j <= i; j++) {
            engine.compute(shells[j], shells[i]); if (result[0] == nullptr) continue;
            Eigen::Map<const Eigen::MatrixXd> buffer(result[0], shells[i].size(), shells[j].size());
            matrix.block(sh2bf[i], sh2bf[j], shells[i].size(), shells[j].size()) = buffer;
            if (i != j) {
                matrix.block(sh2bf[j], sh2bf[i], shells[j].size(), shells[i].size()) = buffer.transpose();
            }
        }
    }
    return matrix;
};

void System::move(Eigen::MatrixXd dir) {
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms.at(i).x += dir(i, 0);
        atoms.at(i).y += dir(i, 1);
        atoms.at(i).z += dir(i, 2);
    }
    shells = libint2::BasisSet(basis, atoms, true);
}
