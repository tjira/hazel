#pragma once

#include "hf.h"

class MP {
public:
    struct OptionsRestricted {
        HF::ResultsRestricted rhfres;
    };
    struct ResultsRestricted {
        double Ecorr;
    };
public:
    // constructor
    MP(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    ResultsRestricted mp2(const System& system, const Tensor<4>& Jmo, bool print = true) const;

private:
    OptionsRestricted ropt;
};
