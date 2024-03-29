#include "transform.h"

Tensor<4> Transform::Coulomb(const Tensor<4>& J, const Matrix& C) {
    // throw an error if coulomb tensor in AO basis is not calculated
    if (!J.size()) throw std::runtime_error("TO TRANSFORM THE COULOMB TENSOR INTO THE MO BASIS YOU NEED TO CALCULATE IT FIRST");

    // declare the ERI tensor in molecular orbital basis and tensors of partial transform
    Tensor<4> J01(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J01.setZero();
    Tensor<4> J02(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J02.setZero();
    Tensor<4> J03(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J03.setZero();
    Tensor<4> Jmo(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); Jmo.setZero();

    // first 5th order partial transform
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < J.dimension(0); i++) {
        for (int a = 0; a < J.dimension(0); a++) {
            for (int b = 0; b < J.dimension(1); b++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J01(i, b, c, d) += J(a, b, c, d) * C(a, i);
                    }
                }
            }
        }
    }

    // second 5th order partial transform
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int b = 0; b < J.dimension(1); b++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J02(i, j, c, d) += J01(i, b, c, d) * C(b, j);
                    }
                }
            }
        }
    }

    // third 5th order partial transform
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int k = 0; k < J.dimension(2); k++) {
                for (int c = 0; c < J.dimension(2); c++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        J03(i, j, k, d) += J02(i, j, c, d) * C(c, k);
                    }
                }
            }
        }
    }

    // fourth 5th order partial transform
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread)
    #endif
    for (int i = 0; i < J.dimension(0); i++) {
        for (int j = 0; j < J.dimension(1); j++) {
            for (int k = 0; k < J.dimension(2); k++) {
                for (int l = 0; l < J.dimension(3); l++) {
                    for (int d = 0; d < J.dimension(3); d++) {
                        Jmo(i, j, k, l) += J03(i, j, k, d) * C(d, l);
                    }
                }
            }
        }
    }

    // return the tensor
    return Jmo;
}

Tensor<4> Transform::CoulombSpin(const Tensor<4>& J, const Matrix& C) {
    // create the spin indices
    Vector ind(2 * C.rows()); std::iota(ind.begin(), ind.end(), 0); ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    BMatrix indm(ind.size(), ind.size()); for (int i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // coefficients in MS and coulomb tensor in AS basis
    Matrix Cs = Utility::Repeat(C, 2, 1).replicate<2, 1>().array() * Utility::Repeat(indm.cast<double>(), C.rows(), 0).topRows(indm.rows()).array();
    Tensor<4> Js = Eigen::Kron<4>(Matrix::Identity(2, 2), Eigen::Kron<4>(Matrix::Identity(2, 2), J).shuffle(Array<4>{3, 2, 1, 0}));

    // return coulomb tensor in MS
    return Coulomb(Js, Cs);
}

Matrix Transform::Oneelec(const Matrix& A, const Matrix& C) {
    // create the transformed matrix
    Matrix Amo(A.rows(), A.cols());

    // perform the transformation
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            for (int a = 0; a < Amo.rows(); a++) {
                for (int b = 0; b < Amo.cols(); b++) {
                    Amo(i, j) += A(a, b) * C(a, i) * C(b, j);
                }
            }
        }
    }

    // return the transformed matrix
    return Amo;
}

Matrix Transform::OneelecSpin(const Matrix& A, const Matrix& C) {
    // expand the dimentsions of the MO basis matrix
    Matrix Ams = Utility::Repeat(Utility::Repeat(Transform::Oneelec(A, C), 2, 0), 2, 1);

    // create the spin indices
    Vector ind(Ams.rows()); std::iota(ind.begin(), ind.end(), 0); ind = ind.unaryExpr([](auto s) {return double(int(s) % 2);});
    BMatrix indm(ind.size(), ind.size()); for (int i = 0; i < ind.size(); i++) indm.row(i) = ind.transpose().array() == ind(i);

    // element wise multiply and return
    return Ams.array() * indm.cast<double>().array();
}
