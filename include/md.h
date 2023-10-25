#pragma once

#include "gradient.h"
#include "hessian.h"

class MD {
public:
    // constructor
    MD(int iters, double step, const std::string& output) : iters(iters), step(step), output(output) {}

    // methods
    void run(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print = true) const;

private:
    int iters; double step; std::string output;
};
