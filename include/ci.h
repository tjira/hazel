#pragma once

#include "data.h"

class CI {
public:
    CI(const Data& data) : data(data) {}

    // order methods
    Data cid(bool print = true) const;

private:
    const Data data;
};
