#include "eigen.h"

template <size_t D> using Ten = Tensor<D>; typedef Matrix Mat;

std::ostream& Eigen::operator<<(std::ostream& os, const Mat& A) {
    os << "(" << A.rows() << "x" << A.cols() << "): ";
    for (int i = 0; i < (A.cols() % 5 ? A.cols() / 5 + 1 : A.cols() / 5); i++) {
        os << i * 5 + 1 << "-" << std::min((int)A.cols(), 5 * i + 5) << "\n";
        for (int j = 0; j < A.rows(); j++) {
            for (int k = 5 * i; k < std::min((int)A.cols(), 5 * i + 5); k++) {
                os << std::setw(20) << A(j, k) << " ";
            }
            if (j < A.rows() - 1) os << "\n";
        }
        if (i < (A.cols() % 5 ? A.cols() / 5 + 1 : A.cols() / 5) - 1) os << "\n";
    }
    return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<3>& A) {
    os << Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1))); return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<4>& A) {
    os << Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3))); return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<5>& A) {
    os << Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2) * A.dimension(4), A.dimension(1) * A.dimension(3))); return os;
}

std::ostream& Eigen::operator<<(std::ostream& os, const Ten<6>& A) {
    os << Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2) * A.dimension(4), A.dimension(1) * A.dimension(3) * A.dimension(5))); return os;
}

std::string Eigen::MemMatrix(const Mat& A) {
    size_t byte = A.rows() * A.cols() * sizeof(double); std::stringstream stream;
    if (byte >= 1e9) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e9 << "GB";
    else if (byte >= 1e6) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e6 << "MB";
    else stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e3 << "KB";
    return stream.str();
}

std::string Eigen::MemTensor(const Tensor<double, 3, Eigen::ColMajor>& A) {
    size_t byte = A.dimension(0) * A.dimension(1) * A.dimension(2) * sizeof(double); std::stringstream stream;
    if (byte >= 1e9) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e9 << "GB";
    else if (byte >= 1e6) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e6 << "MB";
    else stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e3 << "KB";
    return stream.str();
}

std::string Eigen::MemTensor(const Tensor<double, 4, Eigen::ColMajor>& A) {
    size_t byte = A.dimension(0) * A.dimension(1) * A.dimension(2) * A.dimension(3) * sizeof(double); std::stringstream stream;
    if (byte >= 1e9) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e9 << "GB";
    else if (byte >= 1e6) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e6 << "MB";
    else stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e3 << "KB";
    return stream.str();
}

std::string Eigen::MemTensor(const Tensor<double, 5, Eigen::ColMajor>& A) {
    size_t byte = A.dimension(0) * A.dimension(1) * A.dimension(2) * A.dimension(3) * A.dimension(4) * sizeof(double); std::stringstream stream;
    if (byte >= 1e9) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e9 << "GB";
    else if (byte >= 1e6) stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e6 << "MB";
    else stream << std::fixed << std::setprecision(2) << std::setw(6) << std::setfill('0') << byte / 1e3 << "KB";
    return stream.str();
}

void Eigen::Write(const std::string& fname, const Mat& A) {
    // open the output file
    std::ofstream file(fname);

    // write the matrix
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            file << std::fixed << std::setprecision(14) << std::setw(20) << A(i, j) << (j < A.cols() - 1 ? " " : "");
        }
        file << std::endl;
    }
}

void Eigen::Write(const std::string& fname, const Ten<4>& A) {
    Write(fname, Mat(Eigen::Map<const Mat>(A.data(), A.dimension(0) * A.dimension(2), A.dimension(1) * A.dimension(3))));
}
