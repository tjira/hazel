#pragma once

#include <libint2/diis.h>
#include <filesystem>

struct System {
    System(const std::string& fname, const std::string& basis, int charge, int multi);

    // properties of the system
    std::vector<libint2::Atom> atoms;
    int electrons, charge, multi;
    libint2::BasisSet shells;
    std::string basis;
};
