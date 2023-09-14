#pragma once

#include "gradient.h"

struct Lambda {
    // energy and gradient functions
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsUnrestricted& uhfopt, const std::vector<double>& gopt, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D);

    // energy functions
    static std::function<double(System)> EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D);
    static std::function<double(System)> EMP2(const HF::OptionsRestricted& rhfopt, Matrix D);
    static std::function<double(System)> EHF(const HF::OptionsRestricted& rhfopt, Matrix D);
};
