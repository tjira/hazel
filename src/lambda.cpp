#include "lambda.h"

std::function<double(System)> Lambda::ECI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, Matrix D) {
    return [rhfopt, excits, D](System system) {
        // calculate the necessary integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform SCF
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // transform matrces to MS basis
        Matrix Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C);
        Tensor<4> Jms = Transform::CoulombSpin(system.ints.J, rhfres.C);

        // calculate the CI energy and return
        return rhfres.E + CI(rhfres).rci(system, excits, Hms, Jms, false).Ecorr;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        // calculate the necessary integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the SCF and return the energy
        return HF(rhfopt).rscf(system, D, false).E;
    };
}

std::function<double(System)> Lambda::EHF(const HF::OptionsUnrestricted& uhfopt, Matrix D) {
    return [uhfopt, D](System system) {
        // calculate the necessary integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the SCF and return the energy
        return HF(uhfopt).uscf(system, D, false).E;
    };
}

std::function<double(System)> Lambda::EMP2(const HF::OptionsRestricted& rhfopt, Matrix D) {
    return [rhfopt, D](System system) {
        // calculate the necessary integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the SCF
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // transform the J to MO basis
        Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);

        // perform the MP calculation and return the energy
        return rhfres.E + MP(rhfres).rmp2(system, Jmo, false);
    };
}

std::function<double(System)> Lambda::EBAGEL(const Bagel::Options& bagelopt) {
    return [bagelopt](System system) {
        // run the calculation
        auto bagelres = Bagel(system, bagelopt).run();

        // return the energy
        return bagelres.E;
    };
}

std::function<double(System)> Lambda::EORCA(const Orca::Options& orcaopt) {
    return [orcaopt](System system) {
        // run the calculation
        auto orcares = Orca(system, orcaopt).run();

        // return the energy
        return orcares.E;
    };
}

std::function<Vector(System)> Lambda::ESBAGEL(const Bagel::Options& bagelopt) {
    return [bagelopt](System system) {
        // run the calculation
        auto bagelres = Bagel(system, bagelopt).run();

        // return the energy
        return bagelres.excs;
    };
}

std::function<Vector(System)> Lambda::ESORCA(const Orca::Options& orcaopt) {
    return [orcaopt](System system) {
        // run the calculation
        auto orcares = Orca(system, orcaopt).run();

        // return the energy
        return orcares.excs;
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGCI(const HF::OptionsRestricted& rhfopt, const std::vector<int>& excits, double gstep, Matrix D) {
    return [rhfopt, excits, gstep, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the SCF
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the CI gradient
        Matrix G = Gradient(gstep).get(system, Lambda::ECI(rhfopt, excits, rhfres.D), false);

        // transform the Hamiltonian and Coulomb tensor to MS basis
        Matrix Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C);
        Tensor<4> Jms = Transform::CoulombSpin(system.ints.J, rhfres.C);

        // perform the CI calculation and get thecorrelation energy
        double Ecorr = CI(rhfres).rci(system, excits, Hms, Jms, false).Ecorr;

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsRestricted& rhfopt, double gstep, Matrix D) {
    return [rhfopt, gstep, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the HF calculation for the curent geometry
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the numerical or analytical gradient
        if (gstep) return std::tuple{rhfres.E, Gradient(gstep).get(system, Lambda::EHF(rhfopt, rhfres.D), false)};
        else {
            // calculate the integral derivatives
            system.dints.dT = Integral::dKinetic(system), system.dints.dV = Integral::dNuclear(system);
            system.dints.dS = Integral::dOverlap(system), system.dints.dJ = Integral::dCoulomb(system);

            // return the energy and analytical gradient
            return std::tuple{rhfres.E, Matrix(Gradient().get(system, rhfres, false) + Integral::dRepulsion(system))};
        }
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGHF(const HF::OptionsUnrestricted& uhfopt, double gstep, Matrix D) {
    return [uhfopt, gstep, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // perform the HF calculation for the curent geometry
        HF::ResultsUnrestricted uhfres = HF(uhfopt).uscf(system, D, false);

        // calculate the numerical or analytical gradient
        if (gstep) return std::tuple{uhfres.E, Gradient(gstep).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)), false)};
        else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGMP2(const HF::OptionsRestricted& rhfopt, double gstep, Matrix D) {
    return [rhfopt, gstep, D](System& system) {
        // calculate the atomic integrals
        system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
        system.ints.S = Integral::Overlap(system), system.ints.J = Integral::Coulomb(system);

        // recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system, D, false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient(gstep).get(system, Lambda::EMP2(rhfopt, rhfres.D), false);
        Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);
        double Ecorr = MP(rhfres).rmp2(system, Jmo, false);

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGBAGEL(const Bagel::Options& bagelopt) {
    return [bagelopt](System& system) {
        // initialize ORCA object
        Bagel bagel(system, bagelopt);

        // enable the gradient
        bagel.enableGradient();

        // run the calculation
        auto bagelres = bagel.run();

        // calculate the numerical or analytical gradient
        return std::tuple{bagelres.E, bagelres.Gs.at(0)};
    };
}

std::function<std::tuple<double, Matrix>(System&)> Lambda::EGORCA(const Orca::Options& orcaopt, double gstep) {
    return [orcaopt, gstep](System& system) {
        // initialize ORCA object
        Orca orca(system, orcaopt);

        // enable the gradient
        orca.enableGradient(gstep);

        // run the calculation
        auto orcares = orca.run();

        // calculate the numerical or analytical gradient
        return std::tuple{orcares.E, orcares.G};
    };
}

std::function<std::tuple<Vector, std::vector<Matrix>>(System&)> Lambda::ESGSBAGEL(const Bagel::Options& bagelopt, const std::vector<int>& targets) {
    return [bagelopt, targets](System& system) {
        // initialize ORCA object
        Bagel bagel(system, bagelopt);

        // enable the gradient
        bagel.enableGradient(targets);

        // run the calculation
        auto bagelres = bagel.run();

        // calculate the numerical or analytical gradient
        return std::tuple{bagelres.excs, bagelres.Gs};
    };
}
