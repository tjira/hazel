#pragma once

#include "logger.h"
#include "timer.h"
#include <boost/format.hpp>
#include <libint2.hpp>

struct MullikenResult {
    Mat DS;
    Vec q;
};

class System {
public:
    // constructors
    System(std::string filename, std::string basis);
    System(System system, Mat q);

    // getters
    std::vector<libint2::Atom> getAtoms() const { return atoms; };
    libint2::Atom getAtom(size_t i) const { return atoms.at(i); };
    libint2::BasisSet getShells() const { return shells; };
    size_t getSize() const { return atoms.size(); }
    int getElectrons() const { return electrons; }
    double getRepulsion() const;

    // computers
    Mat integralSingle(libint2::Operator op) const;
    Mat integralCoulomb(Mat D) const;
    void move(Mat dir);

private:
    std::vector<libint2::Atom> atoms;
    libint2::BasisSet shells;
    std::string basis;
    int electrons;
};
