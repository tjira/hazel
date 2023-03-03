#include "../include/potential.h"

Eigen::Vector3d AngleHarmonic::F(const Particle& a, const Particle& b, const Particle& c, int i) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords(), q3 = c.getCoords();
    if (i == 1 || i == 3) {
        ret
    }
    return {0, 0, 0};
}

double AngleHarmonic::U(const Particle& a, const Particle& b, const Particle& c) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords(), q3 = c.getCoords();
    double phi = std::acos((q1 - q2).dot((q3 - q2)) / (q1 - q2).norm() / (q3 - q2).norm());
    return k * std::pow((phi - phi0), 2);
}

Eigen::Vector3d BondHarmonic::F(const Particle& a, const Particle& b, int i) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 2 * k * ((q1 - q2).norm() - r0) / (q1 - q2).norm() * (q1 - q2) * (i == 1 ? 1 : -1);
}

double BondHarmonic::U(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return k * std::pow(((q1 - q2).norm() - r0), 2);
}

Eigen::Vector3d LennardJones::F(const Particle& a, const Particle& b, int i) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 24 * epsilon * (2 * std::pow(sigma / (q1 - q2).norm(), 12) - std::pow(sigma / (q1 - q2).norm(), 6)) / std::pow((q1 - q2).norm(), 2) * (q2 - q1) * (i == 1 ? 1 : -1);
}

double LennardJones::U(const Particle& a, const Particle& b) const {
    Eigen::Vector3d q1 = a.getCoords(), q2 = b.getCoords();
    return 4 * epsilon * (std::pow(sigma / (q1 - q2).norm(), 12) - std::pow(sigma / (q1 - q2).norm(), 6));
}
