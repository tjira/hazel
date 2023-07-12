#pragma once

#include "data.h"

class HF {
public:
    // constructor
    HF(const Data& data);

    // methods
    Data scf(bool print = true) const;

private:
    Data data;
};
