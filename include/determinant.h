#pragma once

#include "utility.h"
#include "eigen.h"
#include "timer.h"

class Determinant {
public:
    // constructors
    Determinant(int norb, const std::vector<int>& a, const std::vector<int>& b);
    Determinant(int norb, int nocca, int noccb);

    // operators
    friend std::ostream& operator<<(std::ostream& os, const Determinant& det);

    // excitation generators
    std::vector<Determinant> full() const;

    // integrals and functions
    double hamilton(const Determinant& det, const Matrix& Hms, const Tensor<4>& Jms) const;

private:
    std::vector<int> a, b; int norb;
};
