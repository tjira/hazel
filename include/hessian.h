#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Hessian {
public:
    Hessian(const Data& data) : data(data) {}
    Data frequency(bool print = true) const;
    Data get(bool print = true) const;

private:
    Data get(const std::function<Data(Data)>& efunc, bool print) const;

private:
    Data data;
};
