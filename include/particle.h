#pragma once

#include "ptable.h"
#include <Eigen/Eigen>

class Particle {
    friend class MolecularDynamics;
public:
    Particle(Eigen::Vector3d q, std::string symbol);
    std::string getSymbol() const { return symbol; }
    Eigen::Vector3d getCoords() const { return q; }

private:
    Eigen::Vector3d q, v, a;
    std::string symbol;
    double mass;
};
