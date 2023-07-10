#pragma once

#include "transform.h"
#include "roothaan.h"
#include "data.h"

class MP {
public:
    MP(const Data& data) : Frequency(this), Gradient(this), Hessian(this), Optimizer(this), data(data) {}
    Data mp2(bool print = true) const;
    class Frequency {
    public:
        Frequency(const MP* mp) : mp(mp) {}
        Data mp2(bool print = true) const;
    private:
        const MP* mp;
    } Frequency;
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
    class Hessian {
    public:
        Hessian(const MP* mp) : mp(mp) {}
        Data mp2(bool print = true) const;
    private:
        Data mp2Analytical(bool print = true) const;
        Data mp2Numerical(bool print = true) const;
    private:
        const MP* mp;
    } Hessian;
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
