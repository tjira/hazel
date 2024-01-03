#include "system.h"

System::System(std::ifstream& stream, const std::string& basis, int charge, int multi) : basis(basis), originbasis(basis), electrons(0), charge(charge), multi(multi) {
    // check for the input file existence
    if (!stream.good()) throw std::runtime_error("SYSTEM FILE DOES NOT EXIST");

    // replace the basis placeholders
    std::replace(this->basis.begin(), this->basis.end(), '*', 's'), std::replace(this->basis.begin(), this->basis.end(), '+', 'p');

    // throw an error if impossible combination of charge and multiplicity
    if (std::abs(charge) % 2 == 0 && multi % 2 == 0) {
        throw std::runtime_error("MOLECULE CAN'T HAVE AN EVEN CHARGE AND MULTIPLICITY AT THE SAME TIME.");
    } else if (std::abs(charge) % 2 == 1 && multi % 2 == 1) {
        throw std::runtime_error("MOLECULE CAN'T HAVE AN ODD CHARGE AND MULTIPLICITY AT THE SAME TIME.");
    }

    // assign the atoms and shift the charge
    atoms = libint2::read_dotxyz(stream); electrons -= charge;

    // assign all the other variables
    coords = Matrix(atoms.size(), 3), dists = Matrix(atoms.size(), atoms.size());
    for (const auto& atom : atoms) electrons += atom.atomic_number;
    shells = libint2::BasisSet(this->basis, atoms, true);

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

Determinant System::det() const {
    return Determinant(ints.S.cols(), (electrons + multi - 1) / 2, (electrons - multi + 1) / 2);
}

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

void System::save(const std::string& fname, std::ios::openmode mode) const {
    // open the file, write number of atoms and set output stream flags
    std::ofstream file(fname, mode); file << atoms.size() << "\n" << fname;
    file << std::fixed << std::setprecision(14) << "\n";

    // print all the atom coordinates
    for (int i = 0; i < coords.rows(); i++) {
        file << an2sm.at(atoms.at(i).atomic_number) << " "
             << std::setw(20) << coords(i, 0) << " "
             << std::setw(20) << coords(i, 1) << " "
             << std::setw(20) << coords(i, 2) << std::endl;
    }
}
