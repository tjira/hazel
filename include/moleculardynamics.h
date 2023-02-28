#pragma once

#include "forcefield.h"
#include <boost/format.hpp>
#include <fstream>

class MolecularDynamics {
    typedef MolecularDynamicsOptions Options;
    typedef MolecularDynamicsResult Result;
public:
    MolecularDynamics(const ForceField& field, Options opt);
    Result run(std::vector<Particle> particles, std::string output) const;

private:
    ForceField field;
    Options opt;
};
