#pragma once

#include <Eigen/Eigen>
#include <libint2.hpp>

class Molecule {
public:
    Molecule(std::string filename, std::string basis);
    Eigen::MatrixXd coulomb(Eigen::MatrixXd D) const;
    Eigen::MatrixXd kinetic() const;
    Eigen::MatrixXd nuclear() const;
    Eigen::MatrixXd overlap() const;
    double nuclearRepulsion() const;

private:
    template <int n> Eigen::MatrixXd integral(libint2::Engine engine, Eigen::MatrixXd D = {}) const;
    std::vector<libint2::Atom> atoms; libint2::BasisSet shells;

public:
    int nel() const {
        return std::accumulate(atoms.begin(), atoms.end(), 0, [](int e, const auto& a) { return e + a.atomic_number; });
    }
};

template<> Eigen::MatrixXd Molecule::integral<1>(libint2::Engine engine, Eigen::MatrixXd D) const;
template<> Eigen::MatrixXd Molecule::integral<2>(libint2::Engine engine, Eigen::MatrixXd D) const;
