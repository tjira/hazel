#include "wavetool.h"

Vector Wavetool::Mulliken(const System& system, const Matrix& D) {
    Matrix q = Eigen::VectorXd::Zero(system.atoms.size()), DS = D.cwiseProduct(system.ints.S);
    for (size_t i = 0, j = 0; i < system.shells.size(); i++) {
        for (size_t k = 0; k < system.shells.at(i).size(); k++) {
            q(system.shells.shell2atom(system.atoms).at(i)) -= DS.colwise().sum()(j + k);
        }
        j += system.shells.at(i).size();
    }
    for (size_t i = 0; i < system.atoms.size(); i++) {
        q(i) += system.atoms.at(i).atomic_number;
    }
    return q;
}
