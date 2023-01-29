#include "../include/molecule.h"

namespace libint2 {
    int threadID() {
        #if defined(_OPENMP)
        return omp_get_thread_num();
        #else
        return 0;
        #endif
    }
}

Molecule::Molecule(std::string filename, std::string basis) : filename(filename), setname(basis) {
    std::ifstream file(filename); atoms = libint2::read_dotxyz(file), shells = libint2::BasisSet(basis, atoms, true);
}

int Molecule::getElectronCount() const {
    return std::accumulate(atoms.begin(), atoms.end(), 0, [](int e, const auto& a) { return e + a.atomic_number; });
}

double Molecule::getRepulsion() const {
    auto value = 0.0;
    for (size_t i = 0; i < atoms.size(); i++) {
        for (size_t j = i + 1; j < atoms.size(); j++) {
            double x = atoms.at(i).x - atoms.at(j).x, y = atoms.at(i).y - atoms.at(j).y, z = atoms.at(i).z - atoms.at(j).z;
            value += atoms.at(i).atomic_number * atoms.at(j).atomic_number / Eigen::Vector3d(x, y, z).norm();
        }
    }
    return value;
}

Eigen::MatrixXd Molecule::integral(libint2::Operator op, Eigen::MatrixXd D) const {
    libint2::Engine engine(op, shells.max_nprim(), shells.max_l());
    if (op == libint2::Operator::nuclear) {
        engine.set_params(libint2::make_point_charges(atoms));
    } else if (op == libint2::Operator::coulomb) {
        return integralDouble(engine, D);
    }
    return integralSingle(engine);
}

Eigen::MatrixXd Molecule::integralDouble(libint2::Engine engine2, Eigen::MatrixXd D) const {
    libint2::Engine engine = engine2;
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(shells.nbf(), shells.nbf());
    const auto& result = engine.results(); auto sh2bf = shells.shell2bf();
    for(size_t i = 0; i < shells.size(); i++) {
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

Eigen::MatrixXd Molecule::integralSingle(libint2::Engine engine) const {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(shells.nbf(), shells.nbf());
    std::vector<libint2::Engine> engines(libint2::nthreads, engine);
    #pragma omp parallel for default(none) num_threads(libint2::nthreads) shared(engines, matrix)
    for (size_t i = 0; i < shells.size(); i++) {
        const auto& result = engines.at(libint2::threadID()).results(); auto sh2bf = shells.shell2bf();
        for (size_t j = 0; j <= i; j++) {
            engines.at(libint2::threadID()).compute(shells[j], shells[i]); if (result[0] == nullptr) continue;
            Eigen::Map<const Eigen::MatrixXd> buffer(result[0], shells[i].size(), shells[j].size());
            matrix.block(sh2bf[i], sh2bf[j], shells[i].size(), shells[j].size()) = buffer;
            if (i != j) {
                matrix.block(sh2bf[j], sh2bf[i], shells[j].size(), shells[i].size()) = buffer.transpose();
            }
        }
    }
    return matrix;
};
