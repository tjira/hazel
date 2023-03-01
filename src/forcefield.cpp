#include "../include/forcefield.h"

ForceField::ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, System system) {
    for (int i = 0; i < system.getAtomCount(); i++) {
        for (int j = 0; j < system.getAtomCount(); j++) {
            for (auto entry : pairOpt.config) {
                if (entry.atoms.at(0) == i + 1 && entry.atoms.at(1) == j + 1) {
                    pair[{ i, j }].reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                    pair[{ j, i }].reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                }
            }
            for (auto entry : bondOpt.config) {
                if (entry.atoms.at(0) == i + 1 && entry.atoms.at(1) == j + 1) {
                    bond[{ i, j }].reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1)));
                    bond[{ j, i }].reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1)));
                }
            }
        }
    }
}

Eigen::Vector3d ForceField::F(const std::vector<Particle>& particles, int i) const {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    for (size_t j = 0; j < particles.size(); j++) if(i != j) {
        if (pair.find({ i, j }) != pair.end()) F += pair.at({ i, j })->F(particles.at(i), particles.at(j));
        if (bond.find({ i, j }) != bond.end()) F += bond.at({ i, j })->F(particles.at(i), particles.at(j));
    }
    return F;
}

double ForceField::U(const std::vector<Particle>& particles) const {
    double U = 0;
    for (size_t i = 1; i < particles.size(); i++) for (size_t j = i - 1; j < i; j++) {
        if (pair.find({ i, j }) != pair.end()) U += pair.at({ i, j })->U(particles.at(i), particles.at(j));
        if (bond.find({ i, j }) != bond.end()) U += bond.at({ i, j })->U(particles.at(i), particles.at(j));
    }
    return U;
}
