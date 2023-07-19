#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Gradient {
public:
    struct Options {
        bool numerical; double step;
    };
    struct Results {
        Matrix G;
    };
public:
    // constructor
    Gradient(const Data& data);

    // methods
    Data get(const System& system, bool print = true) const;

private:
    Data get(const System& system, const std::function<Data(System, Data)>& efunc, bool print) const;
    Data getHF(const System& system, bool print) const;

private:
    Data data;
};
