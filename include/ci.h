#pragma once

#include "hf.h"

class CI {
public:
    struct OptionsRestricted {
        HF::ResultsRestricted rhfres;
    };
    struct ResultsRestricted {
        Matrix C, H; Vector eig; double Ecorr;
    };
public:
    CI(const OptionsRestricted& ropt) : ropt(ropt) {}

    // order methods
    ResultsRestricted rcisd(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;
    ResultsRestricted rcid(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;
    ResultsRestricted rcis(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;
    ResultsRestricted rfci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

private:
    // private functions
    ResultsRestricted rsolve(const std::vector<Determinant>& dets, const Matrix& Hms, const Tensor<4>& Jms, bool print = true) const;

    // private variables
    OptionsRestricted ropt;
};
