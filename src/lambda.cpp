#include "lambda.h"

std::function<double(System)> Lambda::EHF(const HF::OptionsRestricted& rhfopt) {
    return [rhfopt](System system) {
        return HF(rhfopt).rscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false).E;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsUnrestricted& uhfopt) {
    return [uhfopt](System system) {
        return HF(uhfopt).uscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false).E;
    };
}

std::function<double(System)> Lambda::EMP2(const HF::OptionsRestricted& rhfopt) {
    return [rhfopt](System system) {
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false);
        return rhfres.E + MP({rhfres}).rmp2(system, Tensor<4>(), false);
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt) {
    return [rhfopt, gopt](System& system) {
        // perform the HF calculation for the curent geometry
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false);

        // calculate the numerical or analytical gradient
        if (gopt.at(0)) return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, Lambda::EHF(rhfopt), false)};
        else return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, rhfres, false)};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt) {
    return [rhfopt, gopt](System& system) {
        // delete the calculated integrals and recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient({gopt.at(1)}).get(system, Lambda::EMP2(rhfopt), false);
        double Ecorr = MP({rhfres}).rmp2(system, Tensor<4>(), false);

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}
