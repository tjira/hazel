#pragma once

#include "molecule.h"
#include "printer.h"
#include "timer.h"
#include <libint2/diis.h>

class HartreeFock {
    typedef HartreeFockOptions Options;
    typedef HartreeFockResult Result;
public:
    HartreeFock(Options opt);
    Result scf(const Molecule& molecule) const;

private:
    Options opt;
};
