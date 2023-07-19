#pragma once

#include "gradient.h"
#include "hessian.h"

class Optimizer {
public:
    struct Options {
        double thresh;
    };
public:
    // constructor
    Optimizer(const Options& opt) : opt(opt) {}

    // methods
    System optimize(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print = true) const;

private:
    Options opt;
};
