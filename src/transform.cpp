#include "transform.h"

Tensor<4> Transform::Coulomb(const Tensor<4>& J, const Matrix& C) {
    // declare the ERI tensor in molecular orbital basis and tensors of partial transform
    Tensor<4> J01(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J01.setZero();
    Tensor<4> J02(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J02.setZero();
    Tensor<4> J03(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); J03.setZero();
    Tensor<4> Jmo(J.dimension(0), J.dimension(1), J.dimension(2), J.dimension(3)); Jmo.setZero();

    // first 5th order partial transform
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) shared(J, J01, C)
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
    #pragma omp parallel for num_threads(nthread) shared(J01, J02, C)
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
    #pragma omp parallel for num_threads(nthread) shared(J02, J03, C)
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
    #pragma omp parallel for num_threads(nthread) shared(J03, Jmo, C)
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

Matrix Transform::Oneelec(const Matrix& A, const Matrix& C) {
    // declare the transformed matrix
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
