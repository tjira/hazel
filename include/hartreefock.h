#pragma once

#include "define.h"
#include "libint.h"
#include "ptable.h"
#include "system.h"

class HartreeFock {
public:
    struct Options {
        struct DIIS {
            int start, keep;
            double damp;
        }; DIIS diis;
        struct GRAD {
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
            bool orben;
        }; PRINT print;
        double thresh; int maxiter;
        bool mulliken;
    };

private:
    struct Flags {
        bool diis = true, silent = false;
    };
    struct GDResult {
        Eigen::MatrixXd G;
        double E;
    };
    struct HFResult {
        Eigen::MatrixXd C, D;
        Eigen::VectorXd eps;
        double E;
    };
    struct MDResult {

    };

public:
    HartreeFock(Options opt) : opt(opt) {};
    MDResult dynamics(System system, Flags flags) const;
    GDResult gradient(System system, Flags flags) const;
    HFResult scf(System system, Flags flags) const;

private:
    const Options opt;
};
