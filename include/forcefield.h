#pragma once

#include "potential.h"

class ForceField {
public:
    ForceField(PotentialArray pair) : pair(pair) {}
    Eigen::Vector3d F(const std::vector<Particle>& particles, int i) const;
    double U(const std::vector<Particle>& particles) const;

private:
    PotentialArray pair;
};
