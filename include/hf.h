#pragma once

#include "integral.h"
#include "system.h"

class HF {
public:
    struct OptionsRestricted {
        struct {int start, keep;} diis;
        double thresh; int maxiter;
        bool nocoulomb;
    };
    struct ResultsRestricted {
        Matrix C, D; Vector eps;
        double E, Eel, Enuc;
    };
    struct OptionsUnrestricted {
        struct {int start, keep;} diis;
        double thresh; int maxiter;
        bool nocoulomb;
    };
    struct ResultsUnrestricted {
        Matrix Ca, Cb, Da, Db;
        double E, Eel, Enuc;
        Vector epsa, epsb;
    };
public:
    // constructor
    HF(const OptionsUnrestricted& uopt) : uopt(uopt) {}
    HF(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    ResultsUnrestricted uscf(const System& system, Matrix D, bool print = true) const;
    ResultsRestricted rscf(const System& system, Matrix D, bool print = true) const;

private:
    OptionsUnrestricted uopt;
    OptionsRestricted ropt;
};
