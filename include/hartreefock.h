#pragma once

#include "ptable.h"
#include "system.h"

class HartreeFock {
public:
    struct Options {
        struct DIIS {
            int start, keep;
            double damp;
            bool on;
        }; DIIS diis;
        struct GRAD {
            struct PRINT {
                bool iter;
            }; PRINT print;
            double increment;
            int nthread;
        }; GRAD engrad;
        struct MDYN {
            std::string output;
            double timestep;
            int steps;
        }; MDYN dyn;
        struct PRINT {
            bool kinetic, oneelec, overlap, density;
            bool orben, mos, iter;
        }; PRINT print;
        double thresh; int maxiter;
        bool mulliken;
    };

private:
    struct GDResult {
        Mat G;
        double E;
    };
    struct HFResult {
        Mat C, D;
        Vec eps;
        double E;
    };
    struct MDResult {

    };

public:
    HartreeFock(Options opt) : opt(opt) {};
    MDResult dynamics(System system, bool silent = false) const;
    GDResult gradient(System system, bool silent = false) const;
    HFResult scf(System system, bool silent = false) const;

private:
    const Options opt;
};
