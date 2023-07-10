#pragma once

#include "data.h"

class Roothaan {
public:
    Roothaan(const Data& data) : data(data) {}
    Data frequency(bool print = true) const;
    Data gradient(bool print = true) const;
    Data optimize(bool print = true) const;
    Data hessian(bool print = true) const;
    Data scf(bool print = true) const;

private:
    Data gradientAnalytical(bool print = true) const;
    Data gradientNumerical(bool print = true) const;
    Data hessianAnalytical(bool print = true) const;
    Data hessianNumerical(bool print = true) const;

private:
    const Data data;
};
