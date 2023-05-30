#pragma once

#include "eigen.h"

struct System {
    System(const std::string& fname, const std::string& basis, int charge, int multi);
    void move(const Matrix& dir);

    // properties of the system
    std::vector<libint2::Atom> atoms;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    Matrix coords, dists;
    std::string basis;
};
