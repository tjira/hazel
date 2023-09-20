#pragma once

#include "gradient.h"
#include "hessian.h"

class QD {
public:
    struct Options {
        std::string potfile;
        int iters, nstates;
        double dt, thresh;
        bool imaginary;
    };
    struct Results {
        std::vector<std::vector<CVector>> states;
        Vector energy; CVector r;
    };
public:
    // constructor
    QD(const Options& opt) : opt(opt) {}

    // methods
    Results run(System system, bool print = true) const;

private:
    Options opt;
};
