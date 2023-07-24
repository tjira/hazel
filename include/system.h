#pragma once

#include "eigen.h"
#include "timer.h"

inline std::unordered_map<int, double> masses = {
    {1, 1.0078400},
    {6, 12.011000},
    {8, 15.999000},
    {9, 18.998403},
};

struct System {
    // constructors and all the functions
    System() = default; System(const std::string& fname, const std::string& basis, int charge, int multi);
    void save(const std::string& fname, std::ios::openmode mode = std::ios::out) const;
    void move(const Matrix& dir); System& clearints();

    // properties of the system
    std::vector<libint2::Atom> atoms;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    Matrix coords, dists;
    std::string basis;

    // containers for integrals
    struct {Tensor<3> dT, dS, dV; Tensor<5> dJ;} dints;
    struct {Matrix T, S, V; Tensor<4> J;} ints;
};
