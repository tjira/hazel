#pragma once

#include "gradient.h"
#include "hessian.h"

class Qdyn {
public:
    struct Options {
        int points, iters, nstates;
        double range, dt, thresh;
        bool imaginary;
    };
    struct Results {
        std::vector<std::vector<CVector>> states;
        Vector energy; CVector r;
    };
public:
    // constructor
    Qdyn(const Options& opt) : opt(opt) {}

    // methods
    Results run(System system, bool print = true) const;

private:
    Options opt;
};
