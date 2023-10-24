#pragma once

#include "hf.h"

class MP {
public:
    // constructor
    MP(const HF::ResultsRestricted& rhfres) : rhfres(rhfres) {}

    // methods
    double rmp2(const System& system, const Tensor<4>& Jmo, bool print = true) const;

private:
    HF::ResultsRestricted rhfres;
};
