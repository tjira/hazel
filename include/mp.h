#pragma once

#include "hf.h"

class MP {
public:
    // constructor
    MP(const Data& data);

    // methods
    Data mp2(bool print = true) const;

private:
    Data data;
};
