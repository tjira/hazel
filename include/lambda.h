#pragma once

#include "bagel.h"
#include "ci.h"
#include "gradient.h"
#include "mp.h"
#include "orca.h"

struct Lambda {
    // energy and gradient functions
    static std::function<std::tuple<double, Matrix>(System&)> EGCI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, double gstep, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsUnrestricted& uhfopt, double gstep, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGMP2(const HF::OptionsRestricted& rhfopt, double gstep, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGHF(const HF::OptionsRestricted& rhfopt, double gstep, Matrix D);
    static std::function<std::tuple<double, Matrix>(System&)> EGBAGEL(const Bagel::Options& bagelopt, double gstep);
    static std::function<std::tuple<double, Matrix>(System&)> EGORCA(const Orca::Options& orcaopt, double gstep);

    // energy functions
    static std::function<double(System)> ECI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, Matrix D);
    static std::function<double(System)> EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D);
    static std::function<double(System)> EMP2(const HF::OptionsRestricted& rhfopt, Matrix D);
    static std::function<double(System)> EHF(const HF::OptionsRestricted& rhfopt, Matrix D);
    static std::function<double(System)> EBAGEL(const Bagel::Options& bagelopt);
    static std::function<double(System)> EORCA(const Orca::Options& orcaopt);

    // multiple energy functions
    static std::function<Vector(System)> ESBAGEL(const Bagel::Options& bagelopt);
    static std::function<Vector(System)> ESORCA(const Orca::Options& orcaopt);
};
