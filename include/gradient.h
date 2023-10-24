#pragma once

#include "hf.h"

class Gradient {
public:
    struct Options {
        double step;
    };
public:
    // constructor
    Gradient(const Options& opt) : opt(opt) {}

    // methods
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print = true) const;
    Matrix get(const System& system, const HF::ResultsRestricted& rhfres, bool print = true) const;

private:
    Options opt;
};
