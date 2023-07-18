#pragma once

#include "data.h"

class HF {
public:
    // constructor
    HF(const Data& data);

    // methods
    Data rscf(bool print = true) const;

private:
    Data data;
};
