#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Hessian {
public:
    struct Options {
        bool numerical; double step;
    };
    struct Results {
        Matrix H; Vector freq;
    };
public:
    // constructor
    Hessian(const Data& data) : data(data) {}

    // methods
    Data frequency(const System& system, bool print = true) const;
    Data get(const System& system, bool print = true) const;

private:
    Data get(const System& system, const std::function<Data(System, Data)>& efunc, bool print) const;

private:
    Data data;
};
