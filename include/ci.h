#pragma once

#include "hf.h"

class CI {
public:
    struct OptionsRestricted {
        HF::ResultsRestricted rhfres;
    };
    struct ResultsRestricted {
        Matrix C, H; Vector eig; double Ecorr;
    };
public:
    CI(const OptionsRestricted& ropt) : ropt(ropt) {}

    // order methods
    ResultsRestricted cid(const System& system, const Tensor<4>& Jmo, bool print = true) const;
    ResultsRestricted cis(const System& system, const Tensor<4>& Jmo, bool print = true) const;

private:
    OptionsRestricted ropt;
};
