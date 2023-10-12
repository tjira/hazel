#pragma once

#include "utility.h"

namespace Transform {
    // transform to MO basis
    Tensor<4> Coulomb(const Tensor<4>& J, const Matrix& C);
    Matrix Oneelec(const Matrix& A, const Matrix& C);

    // transform to MS basis
    Tensor<4> CoulombSpin(const Tensor<4>& J, const Matrix& C);
    Matrix OneelecSpin(const Matrix& A, const Matrix& C);
}
