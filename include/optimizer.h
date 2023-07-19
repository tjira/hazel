#pragma once

#include "gradient.h"
#include "hessian.h"

template <class M>
class Optimizer {
public:
    struct Options {
        double thresh;
    };
    struct Results {
        System system; Matrix G;
    };
public:
    // constructor
    Optimizer(const Data& data) : data(data) {}

    // methods
    Results optimize(const System& system, bool print = true) const;

private:
    Results optimize(const System& system, const std::function<Data(System&, Data)>& egfunc, bool print) const;

private:
    Data data;
};
