#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

typedef Eigen::Vector<double, Eigen::Dynamic> Vector; typedef Eigen::IndexPair<int> Pair;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> Matrix;

template <size_t D> using Tensor = Eigen::Tensor<double, D, Eigen::ColMajor>;
template <size_t D> using Axes = Eigen::array<Eigen::IndexPair<int>, D>;
template <size_t D> using Index = Eigen::array<Eigen::Index, D>;
template <size_t D> using Array = Eigen::array<int, D>;

namespace Eigen {
    std::ostream& operator<<(std::ostream& os, const Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 3, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 4, Eigen::ColMajor>& A);
    std::ostream& operator<<(std::ostream& os, const Tensor<double, 5, Eigen::ColMajor>& A);
}

inline Tensor<2> toTensor(Matrix A) {return Eigen::TensorMap<Tensor<2>>(A.data(), A.rows(), A.cols());}
inline Tensor<1> toTensor(Vector A) {return Eigen::TensorMap<Tensor<1>>(A.data(), A.size());}

inline Matrix toMatrix(Tensor<2> A) {return Eigen::Map<Matrix>(A.data(), A.dimension(0), A.dimension(1));}
inline Vector toVector(Tensor<1> A) {return Eigen::Map<Vector>(A.data(), A.dimension(0));}

#include <iomanip>
