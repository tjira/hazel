#pragma once

#include "potential.h"
#include "particle.h"
#include "forward.h"
#include <fstream>
#include <string>

class MolecularDynamics {
    typedef MolecularDynamicsOptions Options;
    typedef MolecularDynamicsResult Result;
public:
    MolecularDynamics(const Potential& pot, Options opt) : pot(&pot), opt(opt) {};
    Result run(std::vector<Particle> particles, std::string output) const;

private:
    const Potential* pot;
    Options opt;
};
