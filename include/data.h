#pragma once

#include "integral.h"
#include "timer.h"

struct Data {
    struct CI {
        Matrix C, H; Vector eig;
        double Ecorr;
    } ci;
    struct MP {
        double Ecorr;
    } mp;
    struct Roothaan {
        struct {int start, keep;} diis;
        double thresh; int maxiter;
        Matrix C, D; Vector eps;
        double E;
        struct Gradient {
            bool numerical;
            double step;
            Matrix G;
        } grad;
    } roothaan;
    Integrals ints, intsmo;
    System system;
};
