#pragma once

#include "gradient.h"
#include "hessian.h"

class MD {
public:
    struct Options {
        int iters; double step;
        std::string output;
    };
    struct Results {

    };
public:
    // constructor
    MD(const Options& opt) : opt(opt) {}

    // methods
    Results run(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print = true) const;

private:
    Options opt;
};
