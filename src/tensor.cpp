#include "tensor.h"

template <>
Matrix Tensor<4>::contract(Matrix A, const Eigen::array<Eigen::IndexPair<int>, 2>& axes) const {
    Eigen::TensorMap<Eigen::Tensor<double, 2, Eigen::ColMajor>> matrix(A.raw().data(), A.rows(), A.cols());
    Eigen::Tensor<double, 2, Eigen::ColMajor> result = M.contract(matrix, axes);
    return Matrix(Map(result.data(), A.rows(), A.cols()));
}


std::ostream& operator<<(std::ostream& os, const Tensor<4>& A) {
    os << A.M; return os;
}
