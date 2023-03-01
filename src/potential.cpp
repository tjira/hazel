#include "../include/potential.h"
#include <Eigen/src/Core/Matrix.h>

Eigen::Vector3d BondHarmonic::F(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 2 * k * ((q1 - q2).norm() - center) / (q1 - q2).norm() * (q1 - q2);
}

double BondHarmonic::U(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return k * std::pow(((q1 - q2).norm() - center), 2);
}

Eigen::Vector3d LennardJones::F(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 24 * epsilon * (2 * std::pow(sigma / (q1 - q2).norm(), 12) - std::pow(sigma / (q1 - q2).norm(), 6)) / std::pow((q1 - q2).norm(), 2) * (q2 - q1);
}

double LennardJones::U(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 4 * epsilon * (std::pow(sigma / (q1 - q2).norm(), 12) - std::pow(sigma / (q1 - q2).norm(), 6));
}
