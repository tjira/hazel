#pragma once

#include "forward.h"
#include "ptable.h"
#include <Eigen/Eigen>
#include <libint2.hpp>

class Molecule {
public:
    // constructor
    Molecule(std::string filename, std::string basis);

    // getters
    int getElectronCount() const;
    double getRepulsion() const;

    // integral computer
    Eigen::MatrixXd integral(libint2::Operator op, Eigen::MatrixXd D = {}) const;

private:
    // private integral coumputers
    Eigen::MatrixXd integralDouble(libint2::Engine engine, Eigen::MatrixXd D) const;
    Eigen::MatrixXd integralSingle(libint2::Engine engine) const;

    // private variables
    std::vector<libint2::Atom> atoms;
    std::string filename, setname;
    libint2::BasisSet shells;
};
