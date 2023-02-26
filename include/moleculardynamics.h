#pragma once

#include "particle.h"
#include "forward.h"
#include <fstream>
#include <string>

class MolecularDynamics {
    typedef MolecularDynamicsOptions Options;
    typedef MolecularDynamicsResult Result;
public:
    MolecularDynamics(Options opt) : opt(opt) {};
    Result run(std::vector<Particle> particles, std::string output) const;

private:
    Options opt;
};
