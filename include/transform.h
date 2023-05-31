#pragma once

#include "eigen.h"

namespace Transform {
    Tensor<4> Coulomb(const Tensor<4>& J, const Matrix& C);
}
