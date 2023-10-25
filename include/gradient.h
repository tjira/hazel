#pragma once

#include "hf.h"

class Gradient {
public:
    // constructors
    Gradient() {}; Gradient(double step) : step(step) {}

    // methods
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print = true) const;
    Matrix get(const System& system, const HF::ResultsRestricted& rhfres, bool print = true) const;

private:
    double step;
};
