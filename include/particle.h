#pragma once

#include "forward.h"
#include "ptable.h"
#include <Eigen/Eigen>

class Particle {
    friend class MolecularDynamics;
    friend class LennardJones;
public:
    Particle(Eigen::Vector3d q, std::string symbol);
    std::string getSymbol() const { return symbol; }
    int getNumber() const { return sm2an.at(symbol); }
    Eigen::Vector3d getCoords() const { return q; }
    double ekin() const { return 0.5 * mass * v.norm() * v.norm(); }

private:
    Eigen::Vector3d q, v, a;
    std::string symbol;
    double mass;
};
