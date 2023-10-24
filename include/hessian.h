#pragma once

#include "system.h"

class Hessian {
public:
    struct Options {
        double step;
    };
public:
    // constructor
    Hessian(const Options& opt) : opt(opt) {}

    // methods
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print = true) const;
    Vector frequency(const System& system, const Matrix& H, bool print = true) const;

private:
    Options opt;
};
