#include "../include/potential.h"

PotentialArray::PotentialArray(std::string name, std::vector<PotentialCoefficients> entries) : array(119, std::vector<std::shared_ptr<Potential>>(119, nullptr)) {
    for (const PotentialCoefficients& entry : entries) {
        array.at(entry.atoms.at(0)).at(entry.atoms.at(1)).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
        array.at(entry.atoms.at(1)).at(entry.atoms.at(0)).reset(new LennardJones(entry.coefs.at(0), entry.coefs.at(1)));
    }
}

Eigen::Vector3d PotentialArray::F(const Particle& a, const Particle& b) const {
    if (array.at(a.getNumber()).at(b.getNumber()) == nullptr) {
        throw std::runtime_error("Force field not initialized for all atoms.");
    }
    return array.at(a.getNumber()).at(b.getNumber())->F(a, b);
}

double PotentialArray::U(const Particle& a, const Particle& b) const {
    if (array.at(a.getNumber()).at(b.getNumber()) == nullptr) {
        throw std::runtime_error("Force field not initialized for all atoms.");
    }
    return array.at(a.getNumber()).at(b.getNumber())->U(a, b);
}

Eigen::Vector3d LennardJones::F(const Particle& a, const Particle& b) const {
    return 24 * epsilon * (2 * std::pow(sigma / (a.q - b.q).norm(), 12) - std::pow(sigma / (a.q - b.q).norm(), 6)) / std::pow((a.q - b.q).norm(), 2) * (b.q - a.q);
}

double LennardJones::U(const Particle& a, const Particle& b) const {
    return 4 * epsilon * (std::pow(sigma / (a.q - b.q).norm(), 12) - std::pow(sigma / (a.q - b.q).norm(), 6));
}
