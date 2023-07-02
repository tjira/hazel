#pragma once

#include "transform.h"
#include "roothaan.h"
#include "data.h"

class MP {
public:
    MP(const Data& data) : Gradient(this), Optimizer(this), data(data) {}
    Data mp2(bool print = true) const;
    class Gradient {
    public:
        Gradient(const MP* mp) : mp(mp) {}
        Data mp2(bool print = true) const;
    private:
        Data mp2Analytical(bool print = true) const;
        Data mp2Numerical(bool print = true) const;
    private:
        const MP* mp;
    } Gradient;
    class Optimizer {
    public:
        Optimizer(const MP* mp) : mp(mp) {}
        Data mp2(bool print = true) const;
    private:
        const MP* mp;
    } Optimizer;

private:
    const Data data;
};
