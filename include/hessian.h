#pragma once

#include "system.h"

class Hessian {
public:
    // constructor
    Hessian(double step) : step(step) {}

    // methods
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print = true) const;
    Vector frequency(const System& system, const Matrix& H, bool print = true) const;

private:
    double step;
};
