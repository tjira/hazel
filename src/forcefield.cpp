#include "../include/forcefield.h"

ForceField::ForceField(std::string name, std::vector<PotentialCoefficients> entries, System system) {
    pair.resize(system.getAtomCount()); for (auto& row : pair) for (size_t i = 0; i < pair.size(); i++) row.emplace_back(new Potential);
    for (int i = 0; i < system.getAtomCount(); i++) {
        for (int j = 0; j <= i; j++) {
            for (auto entry : entries) {
                bool e12 = entry.atoms.at(0) == system.getAtom(i).atomic_number && entry.atoms.at(1) == system.getAtom(j).atomic_number;
                bool e21 = entry.atoms.at(1) == system.getAtom(i).atomic_number && entry.atoms.at(0) == system.getAtom(j).atomic_number;
                if (e12 || e21) {
                    pair.at(i).at(j).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                    pair.at(j).at(i).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
                }
            }
        }
    }
}

Eigen::Vector3d ForceField::F(const std::vector<Particle>& particles, int i) const {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    for (size_t j = 0; j < particles.size(); j++) if(i != j) {
        F += pair.at(i).at(j)->F(particles.at(i), particles.at(j));
    }
    return F;
}

double ForceField::U(const std::vector<Particle>& particles) const {
    double U = 0;
    for (size_t i = 1; i < particles.size(); i++) for (size_t j = i - 1; j < i; j++) {
        U += pair.at(i).at(j)->U(particles.at(i), particles.at(j));
    }
    return U;
}
