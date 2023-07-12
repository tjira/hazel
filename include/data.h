#pragma once

#include "integral.h"
#include "system.h"
#include "timer.h"

struct Data {
    bool nocoulomb = false;
    struct CI {
        Matrix C, H; Vector eig;
        double Ecorr;
    } ci;
    struct MP {
        double Ecorr;
        struct Frequency {
            bool numerical; double step;
            Vector freq; Matrix H;
        } freq;
        struct Gradient {
            bool numerical;
            double step;
            Matrix G;
        } grad;
        struct Optimizer {
            double thresh;
        } opt;
    } mp;
    struct HF {
        struct {int start, keep;} diis;
        double thresh; int maxiter;
        Matrix C, D; Vector eps;
        double E;
        struct Frequency {
            bool numerical; double step;
            Vector freq; Matrix H;
        } freq;
        struct Gradient {
            bool numerical;
            double step;
            Matrix G;
        } grad;
        struct Optimizer {
            double thresh;
        } opt;
    } hf;
    Integrals ints, intsmo;
    System system;

public:
    Data noints() const {
        Data out = *this; out.ints = {}, out.intsmo = {}; return out;
    }
};

#include "transform.h"
