#pragma once


#include "eigen.h"
#include "system.h"
#include "timer.h"

class Roothaan {
public:
    Roothaan(const System& system, int maxiter, double thresh, const std::pair<int, int>& diis);
    
    // solvers
    Matrix gradient(const Tensor<3>& dT, const Tensor<3>& dV, const Tensor<5>& dJ, const Tensor<3>& dS, const Matrix& C, const Vector& eps) const;
    std::tuple<Matrix, Vector, double> scf(const Matrix& H, const Tensor<4>& J, const Matrix& S, Matrix& D, bool print = true) const;

private:
    struct {int start, keep;} diis;
    double thresh; int maxiter;
    const System system;
};
