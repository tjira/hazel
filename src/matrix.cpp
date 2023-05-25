#include "matrix.h"

std::ostream& operator<<(std::ostream& os, const Matrix& A) {
    for (int i = 0; i < std::max((A.cols() + 1) / 5, 1); i++) {
        if (A.cols() > 5) os << i * 5 + 1 << "-" << std::min(A.cols(), i * 5 + 5) << "\n";
        for (int j = 0; j < A.rows(); j++) {
            for (int k = 5 * i; k < std::min(A.cols(), 5 * i + 5); k++) {
                os << std::format("{:20.14f} ", A(j, k));
            }
            if (j < A.rows() - 1) os << "\n";
        }
        if (i < A.cols() / 5) os << "\n";
    }
    return os;
}

std::tuple<Matrix, Matrix> Matrix::eigensolve(const Matrix& A) const {
    Eigen::GeneralizedSelfAdjointEigenSolver<EigenMat> sol(M, A.M);
    return {Matrix(sol.eigenvectors()), Matrix(sol.eigenvalues())};
}
