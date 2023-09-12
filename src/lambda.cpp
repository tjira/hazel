#include "lambda.h"

std::function<double(System)> Lambda::EHF(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        return HF(rhfopt).rscf(system.clearints(), D, false).E;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D) {
    return [uhfopt, D](System system) {
        return HF(uhfopt).uscf(system.clearints(), D, false).E;
    };
}

std::function<double(System)> Lambda::EMP2(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), D, false);
        return rhfres.E + MP({rhfres}).rmp2(system, Tensor<4>(), false);
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D) {
    return [rhfopt, gopt, D](System& system) {
        // perform the HF calculation for the curent geometry
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), D, false);

        // calculate the numerical or analytical gradient
        if (gopt.at(0)) return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D), false)};
        else return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, rhfres, false)};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D) {
    return [rhfopt, gopt, D](System& system) {
        // delete the calculated integrals and recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), D, false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient({gopt.at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D), false);
        double Ecorr = MP({rhfres}).rmp2(system, Tensor<4>(), false);

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}
