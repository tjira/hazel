#pragma once

#include "determinant.h"

#include <libint2/diis.h>

inline std::unordered_map<int, double> masses = {
    {1, 01.007840},
    {6, 12.011000},
    {7, 14.006700},
    {8, 15.999000},
    {9, 18.998403},
};

inline std::unordered_map<int, std::string> an2sm = {
    { 1,  "H"},
    { 6,  "C"},
    { 7,  "N"},
    { 8,  "O"},
    { 9,  "F"},
    {17, "Cl"}
};

struct System {
    // constructors and all the functions
    System() = default; System(std::ifstream& stream, const std::string& basis, int charge, int multi);
    void save(const std::string& fname, std::ios::openmode mode = std::ios::out) const;
    void move(const Matrix& dir); Determinant det() const;

    // properties of the system
    std::vector<libint2::Atom> atoms;
    std::string basis, originbasis;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    Matrix coords, dists;

    // containers for integrals
    struct {Tensor<3> dT, dS, dV; Tensor<5> dJ;} dints;
    struct {Matrix T, S, V; Tensor<4> J;} ints;
};
