#include "eigen.h"

template <size_t D> using Ten = Tensor<D>; typedef Matrix Mat;

std::ostream& Eigen::operator<<(std::ostream& os, const Mat& A) {
    os << "(" << A.rows() << "x" << A.cols() << "): ";
    for (int i = 0; i < A.cols() / 5 + 1; i++) {
        os << i * 5 + 1 << "-" << std::min((int)A.cols(), 5 * i + 5) << "\n";
        for (int j = 0; j < A.rows(); j++) {
            for (int k = 5 * i; k < std::min((int)A.cols(), 5 * i + 5); k++) {
                os << std::setw(20) << A(j, k) << " ";
            }
            if (j < A.rows() - 1) os << "\n";
        }
        if (i < A.cols() / 5) os << "\n";
    }
    return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<3>& A) {
    os << Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1))); return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<4>& A) {
    int m = A.dimension(0), n = A.dimension(1), o = A.dimension(2), p = A.dimension(3);
    os << Mat(Eigen::Map<const Mat>(A.data(), m * o, n * p)); return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<5>& A) {
    int m = A.dimension(0), n = A.dimension(1), o = A.dimension(2), p = A.dimension(3), q = A.dimension(4);
    os << Mat(Eigen::Map<const Mat>(A.data(), m * o * q, n * p)); return os;
}
