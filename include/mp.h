#pragma once

#include "hf.h"

class MP {
public:
    MP(const Data& data) : data(data) {}
    Data mp2(bool print = true) const;

private:
    const Data data;
};
