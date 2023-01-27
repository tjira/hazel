#pragma once

#include "forward.h"
#include "ptable.h"
#include <Eigen/Eigen>
#include <libint2.hpp>

class Molecule {
public:
    Molecule(std::string filename, std::string basis);
    template <int n>
    Eigen::MatrixXd integral(libint2::Operator op, Eigen::MatrixXd D = {}) const;
    double mass() const; int nel() const;
    double nuclearRepulsion() const;

private:
    libint2::Engine makeEngine(libint2::Operator op) const;
    std::vector<libint2::Atom> atoms;
    std::string filename, setname;
    libint2::BasisSet shells;
};

template<> Eigen::MatrixXd Molecule::integral<1>(libint2::Operator op, Eigen::MatrixXd D) const;
template<> Eigen::MatrixXd Molecule::integral<2>(libint2::Operator op, Eigen::MatrixXd D) const;
