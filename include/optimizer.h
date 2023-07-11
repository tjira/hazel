#pragma once

#include "gradient.h"
#include "hessian.h"

template <class M>
class Optimizer {
public:
    Optimizer(const Data& data) : data(data) {}
    Data optimize(bool print = true) const;

private:
    Data optimize(const std::function<Data(Data)>& egfunc, bool print) const;

private:
    Data data;
};
