#pragma once

#include "constant.h"

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/FFT>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> CMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> BMatrix;
typedef Eigen::Vector<std::complex<double>, Eigen::Dynamic> CVector;
typedef Eigen::Vector<double, Eigen::Dynamic> Vector;

template <size_t D> using CTensor = Eigen::Tensor<std::complex<double>, D, Eigen::ColMajor>;
template <size_t D> using Tensor = Eigen::Tensor<double, D, Eigen::ColMajor>;

template <size_t D> using Axes = Eigen::array<Eigen::IndexPair<int>, D>;
template <size_t D> using Index = Eigen::array<Eigen::Index, D>;
template <size_t D> using Array = Eigen::array<int, D>;

namespace Eigen {
    // ostream operators
    std::ostream& operator<<(std::ostream& os, const Matrix<double, Dynamic, Dynamic, ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 3, ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 4, ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 5, ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 6, ColMajor>& A);

    // math functions
    template <size_t T> Tensor<double, T, ColMajor> Kron(const Matrix<double, Dynamic, Dynamic, ColMajor>& A, const Tensor<double, T, ColMajor>& B);

    // custom functions
    void Write(const std::string& fname, const Matrix<double, Dynamic, Dynamic, ColMajor>& A);
    void Write(const std::string& fname, const Tensor<double, 4, ColMajor>& A);

    // complex functions
    template <typename T> T Conj(const T& A) {return A.unaryExpr([](auto x) {return std::conj(x);});}

    // transforms
    inline CVector fftINV(const CVector& A) {FFT<double> fft; CVector B(A.size()); fft.inv(B, A); return B;}
    inline CVector fftFWD(const CVector& A) {FFT<double> fft; CVector B(A.size()); fft.fwd(B, A); return B;}

    // memory getters
    std::string MemMatrix(const Matrix<double, Dynamic, Dynamic, ColMajor>& A);
    std::string MemTensor(const Tensor<double, 3, ColMajor>& A);
    std::string MemTensor(const Tensor<double, 4, ColMajor>& A);
    std::string MemTensor(const Tensor<double, 5, ColMajor>& A);
}

inline Tensor<2> toTensor(Matrix A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}
inline Tensor<1> toTensor(Vector A) {return Eigen::TensorMap<Tensor<1>>(A.data(), A.size());}

inline Matrix toMatrix(Tensor<4> A) {return Eigen::Map<Matrix>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3));}
inline Matrix toMatrix(Tensor<2> A) {return Eigen::Map<Matrix>(A.data(), A.dimension(0), A.dimension(1));}
inline Vector toVector(Tensor<1> A) {return Eigen::Map<Vector>(A.data(), A.dimension(0));}

template <size_t T>
Tensor<T> Eigen::Kron(const Matrix<double, Dynamic, Dynamic, ColMajor>& A, const Tensor<double, T, ColMajor>& B) {
    Tensor<double, T, ColMajor> C(B.dimension(0), B.dimension(1), A.rows() * B.dimension(2), A.cols() * B.dimension(3));
    for (int i = 0; i < C.dimension(0); i++) {
        for (int j = 0; j < C.dimension(1); j++) {
            for (int k = 0; k < C.dimension(2); k++) {
                for (int l = 0; l < C.dimension(3); l++) {
                    C(i, j, k, l) = A(k / B.dimension(2), l / B.dimension(3)) * B(i, j, k % B.dimension(2), l % B.dimension(3));
                }
            }
        }
    }
    return C;
}

#include <libint2/diis.h>
#include "argparse.hpp"
