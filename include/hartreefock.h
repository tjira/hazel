#pragma once

#include "system.h"
#include "timer.h"
#include <boost/format.hpp>
#include <libint2/diis.h>

class HartreeFock {
    typedef HartreeFockOptions Options;
    typedef HartreeFockResult Result;
public:
    HartreeFock(Options opt);
    Result scf(const System& system) const;

private:
    Options opt;
};
