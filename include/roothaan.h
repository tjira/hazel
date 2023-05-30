#pragma once

#include "integral.h"
#include "timer.h"

class Roothaan {
public:
    Roothaan(const System& system, int maxiter, double thresh, const std::pair<int, int>& diis);
    
    // solvers
    std::tuple<Matrix, Vector, double> scf(const Integrals& ints, Matrix& D, bool print = true) const;
    Matrix gradient(const Integrals& ints, const Matrix& C, const Vector& eps) const;
    std::tuple<System, Integrals, Matrix, Matrix> optimize(Integrals ints) const;

private:
    struct {int start, keep;} diis;
    double thresh; int maxiter;
    const System system;
};
