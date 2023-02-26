#pragma once

#include "system.h"
#include "timer.h"

class HartreeFock {
    typedef HartreeFockOptions Options;
    typedef HartreeFockResult Result;
public:
    HartreeFock(Options opt);
    Result scf(const System& system) const;

private:
    Options opt;
};
