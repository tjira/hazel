#pragma once

#include "hf.h"

class CI {
public:
    struct OptionsRestricted {
        HF::ResultsRestricted rhfres;
        std::vector<int> excits;
    };
    struct ResultsRestricted {
        Matrix C, H; Vector eig; double Ecorr;
    };
public:
    CI(const OptionsRestricted& ropt) : ropt(ropt) {}

    // order methods
    ResultsRestricted rfci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;
    ResultsRestricted rci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

private:
    // private functions
    ResultsRestricted rsolve(const std::vector<Determinant>& dets, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

    // private variables
    OptionsRestricted ropt;
};
