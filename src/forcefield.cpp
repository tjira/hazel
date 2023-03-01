#include "../include/forcefield.h"

ForceField::ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, System system) {
    // create the potential matrices and set all to return zero potentials and forces
    pair.resize(system.getAtomCount()); for (auto& row : pair) for (size_t i = 0; i < pair.size(); i++) row.emplace_back(new Potential);
    bond.resize(system.getAtomCount()); for (auto& row : bond) for (size_t i = 0; i < bond.size(); i++) row.emplace_back(new Potential);

    // assign provided potentials
    for (int i = 0; i < system.getAtomCount(); i++) {
        for (int j = 0; j < system.getAtomCount(); j++) {
            for (auto entry : pairOpt.config) {
                if (entry.atoms.at(0) == i + 1 && entry.atoms.at(1) == j + 1) {
                    pair.at(i).at(j).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                    pair.at(j).at(i).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                }
            }
            for (auto entry : bondOpt.config) {
                if (entry.atoms.at(0) == i + 1 && entry.atoms.at(1) == j + 1) {
                    bond.at(i).at(j).reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1)));
                    bond.at(j).at(i).reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1)));
                }
            }
        }
    }
}

Eigen::Vector3d ForceField::F(const std::vector<Particle>& particles, int i) const {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    for (size_t j = 0; j < particles.size(); j++) if(i != j) {
        F += pair.at(i).at(j)->F(particles.at(i), particles.at(j));
        F += bond.at(i).at(j)->F(particles.at(i), particles.at(j));
    }
    return F;
}

double ForceField::U(const std::vector<Particle>& particles) const {
    double U = 0;
    for (size_t i = 1; i < particles.size(); i++) for (size_t j = i - 1; j < i; j++) {
        U += pair.at(i).at(j)->U(particles.at(i), particles.at(j));
        U += bond.at(i).at(j)->U(particles.at(i), particles.at(j));
    }
    return U;
}
