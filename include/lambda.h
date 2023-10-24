#pragma once

#include "ci.h"
#include "gradient.h"
#include "mp.h"

struct Lambda {
    // energy and gradient functions
    static std::function<std::tuple<double, Matrix>(System&)> EGCI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, const std::vector<double>& gopt, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsUnrestricted& uhfopt, const std::vector<double>& gopt, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D);

    // energy functions
    static std::function<double(System)> ECI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, Matrix D);
    static std::function<double(System)> EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D);
    static std::function<double(System)> EMP2(const HF::OptionsRestricted& rhfopt, Matrix D);
    static std::function<double(System)> EHF(const HF::OptionsRestricted& rhfopt, Matrix D);
};
