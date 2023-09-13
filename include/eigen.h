#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define CFREQ 5140.486777894163
#define BOHR2A 0.529177249
#define A2BOHR 1.889725989

inline int nthread = 1;

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/FFT>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;
typedef Eigen::Vector<double, Eigen::Dynamic> Vector;

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> CMatrix;
typedef Eigen::Vector<std::complex<double>, Eigen::Dynamic> CVector;

template <size_t D> using Tensor = Eigen::Tensor<double, D, Eigen::ColMajor>;
template <size_t D> using Axes = Eigen::array<Eigen::IndexPair<int>, D>;
template <size_t D> using Index = Eigen::array<Eigen::Index, D>;
template <size_t D> using Array = Eigen::array<int, D>;

template <size_t D> using CTensor = Eigen::Tensor<std::complex<double>, D, Eigen::ColMajor>;

namespace Eigen {
    std::ostream& operator<<(std::ostream& os, const Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 3, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 4, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 5, Eigen::ColMajor>& A);

    // custom functions
    void Write(const std::string& fname, const Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A);
    void Write(const std::string& fname, const Tensor<double, 4, Eigen::ColMajor>& A);

    // complex functions
    template <typename T> T Conj(const T& A) {return A.unaryExpr([](auto x) {return std::conj(x);});}

    // transforms
    inline CVector fftINV(const CVector& A) {Eigen::FFT<double> fft; CVector B(A.size()); fft.inv(B, A); return B;}
    inline CVector fftFWD(const CVector& A) {Eigen::FFT<double> fft; CVector B(A.size()); fft.fwd(B, A); return B;}

    // memory getters
    std::string MemMatrix(const Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A);
    std::string MemTensor(const Tensor<double, 3, Eigen::ColMajor>& A);
    std::string MemTensor(const Tensor<double, 4, Eigen::ColMajor>& A);
    std::string MemTensor(const Tensor<double, 5, Eigen::ColMajor>& A);
}

inline Tensor<2> toTensor(Matrix A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}
inline Tensor<1> toTensor(Vector A) {return Eigen::TensorMap<Tensor<1>>(A.data(), A.size());}

inline Matrix toMatrix(Tensor<2> A) {return Eigen::Map<Matrix>(A.data(), A.dimension(0), A.dimension(1));}
inline Vector toVector(Tensor<1> A) {return Eigen::Map<Vector>(A.data(), A.dimension(0));}

#include <libint2/diis.h>
#include "argparse.hpp"
