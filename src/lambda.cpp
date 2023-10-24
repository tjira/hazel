#include "lambda.h"

std::function<double(System)> Lambda::ECI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, Matrix D) {
    return [rhfopt, excits, D](System system) {
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
        Matrix Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C);
        Tensor<4> Jms = Transform::CoulombSpin(system.ints.J, rhfres.C);
        return rhfres.E + CI(rhfres).rci(system, excits, Hms, Jms, false).Ecorr;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);
        return HF(rhfopt).rscf(system, D, false).E;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D) {
    return [uhfopt, D](System system) {
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);
        return HF(uhfopt).uscf(system, D, false).E;
    };
}

std::function<double(System)> Lambda::EMP2(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);
        Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);
        return rhfres.E + MP(rhfres).rmp2(system, Jmo, false);
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGCI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, const std::vector<double>& gopt, Matrix D) {
    return [rhfopt, excits, gopt, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient({gopt.at(1)}).get(system, Lambda::ECI(rhfopt, excits, rhfres.D), false);
        Matrix Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C);
        Tensor<4> Jms = Transform::CoulombSpin(system.ints.J, rhfres.C);
        double Ecorr = CI(rhfres).rci(system, excits, Hms, Jms, false).Ecorr;

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D) {
    return [rhfopt, gopt, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the HF calculation for the curent geometry
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the numerical or analytical gradient
        if (gopt.at(0)) return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D), false)};
        else {
            system.dints.dT = Integral::dKinetic(system), system.dints.dV = Integral::dNuclear(system);
            system.dints.dS = Integral::dOverlap(system), system.dints.dJ = Integral::dCoulomb(system);
            return std::tuple{rhfres.E, Gradient({gopt.at(1)}).get(system, rhfres, false)};
        }
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsUnrestricted& uhfopt, const std::vector<double>& gopt, Matrix D) {
    return [uhfopt, gopt, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the HF calculation for the curent geometry
        HF::ResultsUnrestricted uhfres = HF(uhfopt).uscf(system, D, false);

        // calculate the numerical or analytical gradient
        if (gopt.at(0)) return std::tuple{uhfres.E, Gradient({gopt.at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)), false)};
        else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGMP2(const HF::OptionsRestricted& rhfopt, const std::vector<double>& gopt, Matrix D) {
    return [rhfopt, gopt, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient({gopt.at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D), false);
        Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);
        double Ecorr = MP(rhfres).rmp2(system, Jmo, false);

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}
