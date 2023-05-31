#include "integral.h"

Tensor<4> Integral::Coulomb(const System& system) {
    libint2::Engine engine(libint2::Operator::coulomb, system.shells.max_nprim(), system.shells.max_l(), 0, 1e-12);
    return Double(engine, system);
}

Tensor<5> Integral::dCoulomb(const System& system) {
    libint2::Engine engine(libint2::Operator::coulomb, system.shells.max_nprim(), system.shells.max_l(), 1, 1e-12);
    return dDouble(engine, system);
}

Matrix Integral::Kinetic(const System& system) {
    libint2::Engine engine(libint2::Operator::kinetic, system.shells.max_nprim(), system.shells.max_l(), 0, 1e-12);
    return Single(engine, system);
}

Tensor<3> Integral::dKinetic(const System& system) {
    libint2::Engine engine(libint2::Operator::kinetic, system.shells.max_nprim(), system.shells.max_l(), 1, 1e-12);
    return dSingle(engine, system);
}

Matrix Integral::Nuclear(const System& system) {
    libint2::Engine engine(libint2::Operator::nuclear, system.shells.max_nprim(), system.shells.max_l(), 0, 1e-12);
    engine.set_params(libint2::make_point_charges(system.atoms)); return Single(engine, system);
}

Tensor<3> Integral::dNuclear(const System& system) {
    libint2::Engine engine(libint2::Operator::nuclear, system.shells.max_nprim(), system.shells.max_l(), 1, 1e-12);
    engine.set_params(libint2::make_point_charges(system.atoms)); return dSingle(engine, system);
}

Matrix Integral::Overlap(const System& system) {
    libint2::Engine engine(libint2::Operator::overlap, system.shells.max_nprim(), system.shells.max_l(), 0, 1e-12);
    return Single(engine, system);
}

Tensor<3> Integral::dOverlap(const System& system) {
    libint2::Engine engine(libint2::Operator::overlap, system.shells.max_nprim(), system.shells.max_l(), 1, 1e-12);
    return dSingle(engine, system);
}

Tensor<4> Integral::Double(libint2::Engine& engine, const System& system) {
    // create a result buffer and matrix and a convinient map
    Tensor<4> tensor(system.shells.nbf(), system.shells.nbf(), system.shells.nbf(), system.shells.nbf());
    const auto& result = engine.results(); std::vector<size_t> sh2bf = system.shells.shell2bf();

    // loop over all non-duplicate tensor elements
    for (size_t i = 0; i < system.shells.size(); i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t k = 0; k <= i; k++) {
                for (size_t l = 0; l <= (i == k ? j : k); l++) {
                    // compute the integral over current shells and skip if empty and extract indices of the current basis functions
                    engine.compute(system.shells.at(i), system.shells.at(j), system.shells.at(k), system.shells.at(l)); if (result[0] == nullptr) continue;
                    size_t bfi = sh2bf.at(i), bfj = sh2bf.at(j), bfk = sh2bf.at(k), bfl = sh2bf.at(l);

                    // TODO: the following loop could probably be replaced by assigning using eigen maps

                    // assign the result to the correct position in the tensor using the 8-fold symmetry
                    for (size_t m = 0, q = 0; m < system.shells.at(i).size(); m++) {
                        for (size_t n = 0; n < system.shells.at(j).size(); n++) {
                            for (size_t o = 0; o < system.shells.at(k).size(); o++) {
                                for (size_t p = 0; p < system.shells.at(l).size(); p++, q++) {
                                    tensor(m + bfi, n + bfj, o + bfk, p + bfl) = result.at(0)[q];
                                    tensor(m + bfi, n + bfj, p + bfl, o + bfk) = result.at(0)[q];
                                    tensor(n + bfj, m + bfi, o + bfk, p + bfl) = result.at(0)[q];
                                    tensor(n + bfj, m + bfi, p + bfl, o + bfk) = result.at(0)[q];
                                    tensor(o + bfk, p + bfl, m + bfi, n + bfj) = result.at(0)[q];
                                    tensor(o + bfk, p + bfl, n + bfj, m + bfi) = result.at(0)[q];
                                    tensor(p + bfl, o + bfk, m + bfi, n + bfj) = result.at(0)[q];
                                    tensor(p + bfl, o + bfk, n + bfj, m + bfi) = result.at(0)[q];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // return the resulting tensor
    return tensor;
};

Tensor<5> Integral::dDouble(libint2::Engine& engine, const System& system) {
    // create the integral engine and create shell map
    std::vector<size_t> sh2bf = system.shells.shell2bf();
    const auto& result = engine.results();

    // create the result tensor
    Tensor<5> tensor(system.shells.nbf(), system.shells.nbf(), system.shells.nbf(), system.shells.nbf(), 3); tensor.setZero();

    // loop over all non-duplicate tensor elements
    for (size_t i = 0; i < system.shells.size(); i++) {
        for (size_t j = 0; j < system.shells.size(); j++) {
            for (size_t k = 0; k < system.shells.size(); k++) {
                for (size_t l = 0; l < system.shells.size(); l++) {
                    // compute the integral over current shells and skip if empty and extract indices of the current basis functions
                    engine.compute(system.shells.at(i), system.shells.at(j), system.shells.at(k), system.shells.at(l)); if (result[0] == nullptr) continue;
                    size_t bfi = sh2bf.at(i), bfj = sh2bf.at(j), bfk = sh2bf.at(k), bfl = sh2bf.at(l);

                    // assign the result to the correct position in the tensor
                    for (size_t m = 0, q = 0; m < system.shells.at(i).size(); m++) {
                        for (size_t n = 0; n < system.shells.at(j).size(); n++) {
                            for (size_t o = 0; o < system.shells.at(k).size(); o++) {
                                for (size_t p = 0; p < system.shells.at(l).size(); p++, q++) {
                                    for (int r = 0; r < tensor.dimension(4); r++) {
                                        tensor(m + bfi, n + bfj, o + bfk, p + bfl, r) = result.at(r)[q];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // return the resulting tensor
    return tensor;
};

Matrix Integral::Single(libint2::Engine& engine, const System& system) {
    // create a result buffer and matrix and a convinient map
    Matrix matrix(system.shells.nbf(), system.shells.nbf());
    std::vector<size_t> sh2bf = system.shells.shell2bf();
    const auto& result = engine.results();

    // loop over all matrix elements in the lower triengle of the matrix
    for (size_t i = 0; i < system.shells.size(); i++) {
        for (size_t j = 0; j <= i; j++) {
            // compute the integral over corresponding shells and skip if zero
            engine.compute(system.shells.at(j), system.shells.at(i)); if (result.at(0) == nullptr) continue;

            // extract the integral result to a buffer
            Eigen::Map<const Matrix> buffer(result.at(0), system.shells.at(i).size(), system.shells.at(j).size());

            // assign the buffer to the correct position in the result matrix
            matrix.block(sh2bf.at(i), sh2bf.at(j), system.shells.at(i).size(), system.shells.at(j).size()) = buffer;

            // if the element is not on the diagonal fill also the mirrored element
            if (i != j) {
                matrix.block(sh2bf.at(j), sh2bf.at(i), system.shells.at(j).size(), system.shells.at(i).size()) = buffer.transpose();
            }
        }
    }

    // return the resulting integral matrix
    return matrix;
};

Tensor<3> Integral::dSingle(libint2::Engine& engine, const System& system) {
    // create result buffer and shell map
    std::vector<size_t> sh2bf = system.shells.shell2bf();
    const auto& result = engine.results();

    // initiaize the result matrix
    int dim = engine.oper() == libint2::Operator::nuclear ? 3 * system.atoms.size() + 6 : 6;
    Tensor<3> tensor(system.shells.nbf(), system.shells.nbf(), dim); tensor.setZero();

    // loop over all matrix elements
    for (size_t i = 0; i < system.shells.size(); i++) {
        for (size_t j = 0; j < system.shells.size(); j++) {
            // compute the integral over corresponding shells and skip if zero
            engine.compute(system.shells.at(i), system.shells.at(j)); if (result.at(0) == nullptr) continue;

            // assign the result to the correct position in the tensor
            for (size_t m = 0, q = 0; m < system.shells.at(i).size(); m++) {
                for (size_t n = 0; n < system.shells.at(j).size(); n++, q++) {
                    for (int o = 0; o < tensor.dimension(2); o++) {
                        tensor(m + sh2bf.at(i), n + sh2bf.at(j), o) = result.at(o)[q];
                    }
                }
            }
        }
    }

    // return the resulting integral tensor
    return tensor;
};

double Integral::Repulsion(const System& system) {
    // create nuclear repulsion variable
    double repulsion = 0;

    // loop over every nuclear pair exactly once
    for (size_t i = 0; i < system.atoms.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            // extract the atoms and their distance
            libint2::Atom a = system.atoms.at(i), b = system.atoms.at(j);
            double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;
            double dist = std::sqrt(x * x + y * y + z * z);

            // add the value to the repulsion
            repulsion += a.atomic_number * b.atomic_number / dist;
        }
    }

    // return the repulsion
    return repulsion;
}

Matrix Integral::dRepulsion(const System& system) {
    // create nuclear repulsion variable
    Matrix repulsion(system.atoms.size(), 3);

    // loop over every nuclear pair exactly once
    for (size_t i = 0; i < system.atoms.size(); i++) {
        for (size_t j = 0; j < system.atoms.size(); j++) {
            // skip if the same atom
            if (i == j) continue;

            // extract the atoms and their distance
            libint2::Atom a = system.atoms.at(i), b = system.atoms.at(j);
            double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;
            double dist = std::sqrt(x * x + y * y + z * z);

            // add the value to the repulsion
            repulsion(i, 0) -= x * a.atomic_number * b.atomic_number / std::pow(dist, 3);
            repulsion(i, 1) -= y * a.atomic_number * b.atomic_number / std::pow(dist, 3);
            repulsion(i, 2) -= z * a.atomic_number * b.atomic_number / std::pow(dist, 3);
        }
    }

    // return the repulsion
    return repulsion;
}
