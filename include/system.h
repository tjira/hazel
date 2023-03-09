#pragma once

#include "logger.h"
#include "timer.h"
#include <boost/format.hpp>
#include <libint2.hpp>

struct MullikenResult {
    Eigen::MatrixXd DS;
    Eigen::VectorXd q;
};

class System {
public:
    // constructor
    System(std::string filename, std::string basis);

    // getters
    std::vector<libint2::Atom> getAtoms() const { return atoms; };

    // computers
    Eigen::MatrixXd integralSingle(libint2::Operator op) const;
    Eigen::MatrixXd integralCoulomb(Eigen::MatrixXd D) const;
    MullikenResult mulliken(Eigen::MatrixXd D) const;
    void move(Eigen::MatrixXd dir);

private:
    std::vector<libint2::Atom> atoms;
    libint2::BasisSet shells;
    std::string basis;
};
