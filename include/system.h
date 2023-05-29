#pragma once

#include <libint2/diis.h>
#include <filesystem>

struct System {
    System(const std::string& fname, const std::string& basis, int charge, int multi);

    // properties of the system
    int electrons, charge, multi, nocc;
    std::vector<libint2::Atom> atoms;
    libint2::BasisSet shells;
    std::string basis;
};
