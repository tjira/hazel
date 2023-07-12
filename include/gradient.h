#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Gradient {
public:
    // constructor
    Gradient(const Data& data);

    // methods
    Data get(bool print = true) const;

private:
    Data get(const std::function<Data(Data)>& efunc, bool print) const;
    Data getHF(bool print) const;

private:
    Data data;
};
