#include "../include/forcefield.h"

Eigen::Vector3d ForceField::F(const std::vector<Particle>& particles, int i) const {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    for (size_t j = 0; j < particles.size(); j++) {
        if(i != j) F += pair.F(particles.at(i), particles.at(j));
    }
    return F;
}

double ForceField::U(const std::vector<Particle>& particles) const {
    double U = 0;
    for (size_t j = 1; j < particles.size(); j++) {
        for (size_t k = j - 1; k < j; k++) {
            U += pair.U(particles.at(j), particles.at(k));
        }
    }
    return U;
}
