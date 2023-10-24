#pragma once

#include "hf.h"

class CI {
public:
    struct ResultsRestricted {
        Matrix C, H; Vector eig; double Ecorr;
    };
public:
    CI(const HF::ResultsRestricted& rhfres) : rhfres(rhfres) {}

    // order methods
    ResultsRestricted rci(const System& system, const std::vector<int>& excits, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;
    ResultsRestricted rfci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

private:
    // private functions
    ResultsRestricted rsolve(const std::vector<Determinant>& dets, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

    // private variables
    HF::ResultsRestricted rhfres;
};
