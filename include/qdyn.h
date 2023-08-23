#pragma once

#include "gradient.h"
#include "hessian.h"

class Qdyn {
public:
    struct Options {
        int iters; double step;
    };
    struct Results {

    };
public:
    // constructor
    Qdyn(const Options& opt) : opt(opt) {}

    // methods
    Results run(System system, bool print = true) const;

private:
    Options opt;
};
