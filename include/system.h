#pragma once

#include "eigen.h"

inline std::unordered_map<int, double> masses = {
    {1, 1.0078400},
    {6, 12.011000},
    {8, 15.999000},
    {9, 18.998403},
};

struct System {
    System() = default; System(const std::string& fname, const std::string& basis, int charge, int multi);
    void move(const Matrix& dir); void save(const std::string& fname) const;

    // properties of the system
    std::vector<libint2::Atom> atoms;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    Matrix coords, dists;
    std::string basis;
};
