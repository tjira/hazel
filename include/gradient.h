#pragma once

#include "ci.h"
#include "hf.h"
#include "mp.h"

template <class M>
class Gradient {
public:
    struct OptionsRestricted {
        HF::OptionsRestricted rhfopt;
        HF::ResultsRestricted rhfres;
        bool numerical; double step;
    };
public:
    // constructor
    Gradient(const OptionsRestricted& ropt) : ropt(ropt) {}

    // methods
    Matrix get(const System& system, bool print = true) const;

private:
    Matrix get(const System& system, const std::function<double(System)>& efunc, bool print) const;
    Matrix getHF(const System& system, bool print) const;

private:
    OptionsRestricted ropt;
};
