#pragma once

#include "integral.h"
#include "system.h"

class HF {
public:
    struct OptionsRestricted {
        struct {int start, keep;} diis;
        double thresh; int maxiter;
    };
    struct ResultsRestricted {
        Matrix C, D; Vector eps;
        double E, Eel, Enuc;
    };
public:
    // constructor
    HF(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    ResultsRestricted rscf(const System& system, Matrix D, bool print = true) const;

private:
    OptionsRestricted ropt;
};
