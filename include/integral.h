#pragma once

#include "matrix.h"
#include "system.h"
#include "tensor.h"

namespace Integral {
    // specific integrals
    Tensor<4> Coulomb(const System& system);
    Matrix Kinetic(const System& system);
    Matrix Overlap(const System& system);
    Matrix Nuclear(const System& system);

    // general integrals
    Tensor<4> Double(libint2::Engine& engine, const System& system);
    Matrix Single(libint2::Engine& engine, const System& system);

    // additional calculations
    double Repulsion(const System& system);
}
