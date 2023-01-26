#pragma once

#include <Eigen/Eigen>
#include <libint2.hpp>

class Molecule {
public:
    Molecule(std::string filename, std::string basis);
    template <int n>
    Eigen::MatrixXd integral(libint2::Operator op, Eigen::MatrixXd D = {}) const;
    double nuclearRepulsion() const;

private:
    libint2::Engine makeEngine(libint2::Operator op) const;
    std::vector<libint2::Atom> atoms;
    std::string filename, setname;
    libint2::BasisSet shells;

public:
    std::string basis() const {
        return setname;
    }
    std::string fname() const {
        return filename;
    }
    int nel() const {
        return std::accumulate(atoms.begin(), atoms.end(), 0, [](int e, const auto& a) { return e + a.atomic_number; });
    }
};

template<> Eigen::MatrixXd Molecule::integral<1>(libint2::Operator op, Eigen::MatrixXd D) const;
template<> Eigen::MatrixXd Molecule::integral<2>(libint2::Operator op, Eigen::MatrixXd D) const;
