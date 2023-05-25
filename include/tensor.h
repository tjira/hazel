#pragma once

#include "matrix.h"

template <int D>
class Tensor {
public:
    Tensor(auto ...args) : M(args...) {}

    // arithmetic operators
    Tensor operator+(const Tensor& A) const {return Tensor(M + A.M);}
    Tensor operator-(const Tensor& A) const {return Tensor(M - A.M);}
    Tensor operator*(double a) const {return Tensor(a * M);}

    // access operators
    double operator()(auto ...args) const {return M(args...);}
    double& operator()(auto ...args) {return M(args...);}

    // other operators
    friend std::ostream& operator<<(std::ostream& os, const Tensor<4>& A);

    // math functions
    Matrix contract(Matrix A, const Eigen::array<Eigen::IndexPair<int>, 2>& axes) const;
    Tensor t(const Eigen::array<int, D>& dim) const {return Tensor(M.shuffle(dim));}

private:
    Eigen::Tensor<double, D, Eigen::ColMajor> M;
};

template <int D>
Tensor<D> operator*(double a, const Tensor<D>& A) {return A * a;}
