#pragma once

#include "matrix.h"
#include "system.h"
#include "tensor.h"
#include "timer.h"

class Roothaan {
public:
    Roothaan(int maxiter, double thresh);
    
    // solvers
    std::tuple<Matrix, Matrix, double> scf(const Matrix& H, const Tensor<4>& J, const Matrix& S, Matrix D, int nocc, bool print = true) const;

private:
    double thresh; int maxiter;
};
