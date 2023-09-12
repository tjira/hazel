#pragma once

#include "gradient.h"

struct Lambda {
    // energy and gradient functions
    static std::function<std::tuple<double, Matrix>(System&)> EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt);

    // energy functions
    static std::function<double(System)> EHF(const HF::OptionsUnrestricted& uhfopt);
    static std::function<double(System)> EMP2(const HF::OptionsRestricted& rhfopt);
    static std::function<double(System)> EHF(const HF::OptionsRestricted& rhfopt);
};
