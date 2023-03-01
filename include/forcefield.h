#pragma once

#include "potential.h"
#include "system.h"

class ForceField {
public:
    ForceField(PotentialOptions pairOpt, PotentialOptions bondOpt, System system);
    Eigen::Vector3d F(const std::vector<Particle>& particles, int i) const;
    double U(const std::vector<Particle>& particles) const;

private:
    std::vector<std::vector<std::shared_ptr<Potential>>> pair, bond;
};
