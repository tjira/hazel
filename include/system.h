#pragma once

#include "forward.h"
#include "particle.h"
#include "ptable.h"
#include <Eigen/Eigen>
#include <libint2.hpp>

class System {
public:
    // constructor
    System(std::string filename, std::string basis = "6-31G");

    // getters
    std::vector<Particle> getParticles() const;
    libint2::Atom getAtom(int i) const;
    double getNuclearRepulsion() const;
    int getElectronCount() const;
    int getAtomCount() const;

    // computers
    Eigen::MatrixXd integral(libint2::Operator op, Eigen::MatrixXd D = {}) const;
    MullikenResult mulliken(Eigen::MatrixXd D) const;

private:
    // private integral coumputers
    Eigen::MatrixXd integralDouble(libint2::Engine engine, Eigen::MatrixXd D) const;
    Eigen::MatrixXd integralSingle(libint2::Engine engine) const;

    // private variables
    std::vector<libint2::Atom> atoms;
    std::string filename, setname;
    libint2::BasisSet shells;
};
