#pragma once

#include "lambda.h"

class Optimizer {
public:
    // constructor
    Optimizer(double thresh) : thresh(thresh) {}

    // methods
    System optimize(System system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print = true) const;

private:
    double thresh;
};
