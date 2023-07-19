#pragma once

#include "data.h"

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
    HF(const Data& data);

    // methods
    Data rscf(const System& system, bool print = true) const;

private:
    Data data;
};
