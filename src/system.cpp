#include "system.h"

System::System(const std::string& fname, const std::string& basis, int charge, int multi) : electrons(0), charge(charge), multi(multi), basis(basis) {
    if (!std::filesystem::exists(fname)) throw std::runtime_error("System file does not exist.");
    std::ifstream file(fname); atoms = libint2::read_dotxyz(file); electrons -= charge;
    for (const auto& atom : atoms) electrons += atom.atomic_number;
    shells = libint2::BasisSet(basis, atoms, true);
    nocc = electrons / 2;
}
