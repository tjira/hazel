#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Hessian {
public:
    struct OptionsRestricted {
        HF::OptionsRestricted rhfopt;
        HF::ResultsRestricted rhfres;
        bool numerical; double step;
    };
public:
    // constructor
    Hessian(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    std::tuple<Matrix, Vector> frequency(const System& system, bool print = true) const;

private:
    // general function
    Matrix get(const System& system, bool print = true) const;

    // specialized functions
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print) const;

private:
    OptionsRestricted ropt;
};
