#pragma once

#include "data.h"

class HF {
public:
    HF(const Data& data) : data(data) {}
    Data scf(bool print = true) const;

private:
    const Data data;
};
