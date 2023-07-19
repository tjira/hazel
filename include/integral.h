#pragma once

#include "transform.h"
#include "system.h"

namespace Integral {
    // specific integrals
    Matrix Coulomb(const System& system, const Matrix& D);
    Tensor<4> Coulomb(const System& system);
    Matrix Kinetic(const System& system);
    Matrix Overlap(const System& system);
    Matrix Nuclear(const System& system);

    // specific first derivative integrals
    Tensor<5> dCoulomb(const System& system);
    Tensor<3> dKinetic(const System& system);
    Tensor<3> dOverlap(const System& system);
    Tensor<3> dNuclear(const System& system);

    // general integrals
    Matrix Double(libint2::Engine& engine, const System& system, const Matrix& D);
    Tensor<4> Double(libint2::Engine& engine, const System& system);
    Matrix Single(libint2::Engine& engine, const System& system);

    // general derivative integrals
    Tensor<5> dDouble(libint2::Engine& engine, const System& system);
    Tensor<3> dSingle(libint2::Engine& engine, const System& system);

    // additional calculations
    Matrix dRepulsion(const System& system);
    double Repulsion(const System& system);
}
