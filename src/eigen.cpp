#include "eigen.h"

template <size_t D> using Ten = Tensor<D>; typedef Matrix Mat;

std::ostream& Eigen::operator<<(std::ostream& os, const Mat& A) {
    for (int i = 0; i < A.cols() / 6 + 1; i++) {
        os << i * 5 + 1 << "-" << std::min((int)A.cols(), 5 * i + 5) << "\n";
        for (int j = 0; j < A.rows(); j++) {
            for (int k = 5 * i; k < std::min((int)A.cols(), 5 * i + 5); k++) {
                os << std::format("{:20.14f} ", A(j, k));
            }
            if (j < A.rows() - 1) os << "\n";
        }
        if (i < A.cols() / 6) os << "\n";
    }
    return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<3>& A) {
    for (int i = 0; i < A.dimension(2); i++) {
        // extract the dimensions
        auto dim = A.dimensions();

        // slice the provided tensor
        Eigen::array<Eigen::Index, 3> start = {0, 0, i}, size = {dim.at(0), dim.at(1), 1};
        Eigen::Tensor<double, 3, Eigen::ColMajor> tensor = A.slice(start, size);

        // convert the slice to a matrix and print it
        Mat matrix = Eigen::Map<Mat>(tensor.data(), dim.at(0), dim.at(1));
        os << i + 1 << ": " << matrix << (i < dim.at(2) - 1 ? "\n" : "");
    }
    return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<4>& A) {
    for (int i = 0; i < A.dimension(3); i++) {
        // extract the dimensions
        auto dim = A.dimensions();

        // slice the provided tensor
        Eigen::array<Eigen::Index, 4> start = {0, 0, 0, i}, size = {dim.at(0), dim.at(1), dim.at(2), 1};
        Eigen::Tensor<double, 4, Eigen::ColMajor> tensor = A.slice(start, size);

        // convert the slice to a 3rd order tensor and print it
        Ten<3> result = Eigen::TensorMap<Ten<3>>(tensor.data(), dim.at(0), dim.at(1), dim.at(2));
        os << i + 1 << ": " << result << (i < dim.at(3) - 1 ? "\n" : "");
    }
    return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<5>& A) {
    for (int i = 0; i < A.dimension(4); i++) {
        // extract the dimensions
        auto dim = A.dimensions();

        // slice the provided tensor
        Eigen::array<Eigen::Index, 5> start = {0, 0, 0, 0, i}, size = {dim.at(0), dim.at(1), dim.at(2), dim.at(3), 1};
        Eigen::Tensor<double, 5, Eigen::ColMajor> tensor = A.slice(start, size);

        // convert the slice to a 3rd order tensor and print it
        Ten<4> result = Eigen::TensorMap<Ten<4>>(tensor.data(), dim.at(0), dim.at(1), dim.at(2), dim.at(3));
        os << i + 1 << ": " << result << (i < dim.at(4) - 1 ? "\n" : "");
    }
    return os;
}
