#include "../include/potential.h"

Eigen::Vector3d LennardJones::F(const Particle& a, const Particle& b) const {
    return 24 * epsilon * (2 * std::pow(sigma / (a.q - b.q).norm(), 12) - std::pow(sigma / (a.q - b.q).norm(), 6)) / std::pow((a.q - b.q).norm(), 2) * (b.q - a.q);
}

double LennardJones::U(const Particle& a, const Particle& b) const {
    return 4 * epsilon * (std::pow(sigma / (a.q - b.q).norm(), 12) - std::pow(sigma / (a.q - b.q).norm(), 6));
}
