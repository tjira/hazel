#pragma once

#include "system.h"


class MP {
public:
    MP(const System& system);

    // order methods
    double mp2(const Tensor<4>& JMO, const Vector& eps) const;

private:
    int nocc;
};
