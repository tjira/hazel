#include "../include/forcefield.h"

ForceField::ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, PotentialOptions angleOpt) {
    for (auto entry : pairOpt.config) {
        pair[{ entry.atoms.at(0) - 1, entry.atoms.at(1) - 1 }].reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
    }
    for (auto entry : bondOpt.config) {
        bond[{ entry.atoms.at(0) - 1, entry.atoms.at(1) - 1 }].reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1)));
    }
    for (auto entry : angleOpt.config) {
        /* angle[{ entry.atoms.at(0) - 1, entry.atoms.at(1) - 1, entry.atoms.at(2) - 1 }].reset(new BondHarmonic(entry.coefs.at(0), entry.coefs.at(1))); */
    }
}

Eigen::Vector3d ForceField::F(const std::vector<Particle>& particles, int i) const {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    for (const auto& [atoms, pot] : pair) { auto [a, b] = atoms;
        if (int c = a == i ? 1 : b == i ? 2 : 0; c) F += pot->F(particles.at(a), particles.at(b), c);
    }
    for (const auto& [atoms, pot] : bond) { auto [a, b] = atoms;
        if (int c = a == i ? 1 : b == i ? 2 : 0; c) F += pot->F(particles.at(a), particles.at(b), c);
    }
    return F;
}

double ForceField::U(const std::vector<Particle>& particles) const {
    double U = 0;
    for (const auto& [atoms, pot] : pair) { auto [a, b] = atoms;
        U += pot->U(particles.at(a), particles.at(b));
    }
    for (const auto& [atoms, pot] : bond) { auto [a, b] = atoms;
        U += pot->U(particles.at(a), particles.at(b));
    }
    for (const auto& [atoms, pot] : angle) { auto [a, b, c] = atoms;
        U += pot->U(particles.at(a), particles.at(b), particles.at(c));
    }
    return U;
}
