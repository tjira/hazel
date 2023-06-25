#include "system.h"

std::unordered_map<int, std::string> an2sm = {
    {1, "H"},
    {6, "C"},
    {8, "O"}
};

System::System(const std::string& fname, const std::string& basis, int charge, int multi) : electrons(0), charge(charge), multi(multi), basis(basis) {
    if (!std::ifstream(fname).good()) throw std::runtime_error("System file does not exist.");
    std::ifstream file(fname); atoms = libint2::read_dotxyz(file); electrons -= charge;
    coords = Matrix(atoms.size(), 3), dists = Matrix(atoms.size(), atoms.size());
    for (const auto& atom : atoms) electrons += atom.atomic_number;
    shells = libint2::BasisSet(basis, atoms, true);

    // fill the coordinate matrix
    for (int i = 0; i < coords.rows(); i++) {
        coords.row(i) = BOHR2A * Eigen::Vector<double, 3>(atoms.at(i).x, atoms.at(i).y, atoms.at(i).z);
    }

    // fill the distance matrix
    for (int i = 0; i < coords.rows(); i++) {
        for (int j = 0; j < coords.rows(); j++) {
            dists(i, j) = (coords.row(i) - coords.row(j)).norm();
        }
    }
}

// Moves the system in the desired direction. The 'dir' variable needs to be in Bohrs.
void System::move(const Matrix& dir) {
    // shift the atoms
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms.at(i).x += dir(i, 0);
        atoms.at(i).y += dir(i, 1);
        atoms.at(i).z += dir(i, 2);
    }

    // shift the coords and recreate the basis
    shells = libint2::BasisSet(basis, atoms, true), coords += BOHR2A * dir;

    // recalculate the distance matrix
    for (int i = 0; i < coords.rows(); i++) {
        for (int j = 0; j < coords.rows(); j++) {
            dists(i, j) = (coords.row(i) - coords.row(j)).norm();
        }
    }
}

void System::save(const std::string& fname) const {
    std::ofstream file(fname); file << fname << "\n" << atoms.size();
    file << std::fixed << std::setprecision(14) << "\n";
    for (int i = 0; i < coords.rows(); i++) {
        file << an2sm.at(atoms.at(i).atomic_number) << " "
              << std::setw(20) << coords(i, 0) << " "
              << std::setw(20) << coords(i, 1) << " "
              << std::setw(20) << coords(i, 2) << " "
              << std::endl;
    }
}
