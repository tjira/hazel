#pragma once

#include "gradient.h"
#include "hessian.h"

template <class M>
class Optimizer {
public:
    struct OptionsRestricted {
        HF::OptionsRestricted rhfopt;
        Gradient<M>::OptionsRestricted gradopt;
        double thresh;
    };
public:
    // constructor
    Optimizer(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    System optimize(const System& system, bool print = true) const;

private:
    System optimize(const System& system, const std::function<std::tuple<double, Matrix>(System&)>& egfunc, bool print) const;

private:
    OptionsRestricted ropt;
};
