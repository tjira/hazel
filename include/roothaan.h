#pragma once

#include "data.h"

class Roothaan {
public:
    Roothaan(const Data& data) : data(data) {}
    Data gradient(bool print = true) const;
    Data optimize(bool print = true) const;
    Data scf(bool print = true) const;

private:
    Data gradientAnalytical(bool print = true) const;
    Data gradientNumerical(bool print = true) const;

private:
    const Data data;
};
