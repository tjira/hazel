#pragma once

#include "integral.h"
#include "system.h"
#include "transform.h"


class CI {
public:
    CI(const System& system, double Eel);

    // order methods
    std::tuple<Matrix, Matrix, Vector> cid(const Integrals& ints, const Tensor<4>& JMO, const Matrix& C) const;

private:
    double Eel; int nocc, nbf;
};
