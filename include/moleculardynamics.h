#pragma once

#include "particle.h"
#include <boost/format.hpp>
#include <fstream>

struct Field {
    Eigen::Vector3d F(const std::vector<Particle>& particles, int i) const { return { 0, 0, 0 }; }
    double U(const std::vector<Particle>& particles) const { return 0; }
};

class MolecularDynamics {
    typedef MolecularDynamicsOptions Options;
    typedef MolecularDynamicsResult Result;
public:
    MolecularDynamics(Options opt);
    Result run(std::vector<Particle> particles, std::string output) const;

private:
    Field field;
    Options opt;
};
