#pragma once

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> EigenMat;
typedef Eigen::Block<const EigenMat> ConstBlock; typedef Eigen::Block<EigenMat> Block;
typedef Eigen::Map<const EigenMat> ConstMap; typedef Eigen::Map<EigenMat> Map;
typedef Eigen::IndexPair<int> Ind;

#define WIDTH 104
#include <format>

class Matrix {
public:
    Matrix(int rows, int cols) : M(EigenMat::Zero(rows, cols)) {}
    Matrix(int size) : M(EigenMat::Zero(size, size)) {}
    Matrix(const EigenMat& A) : M(A) {};

    // arithmetic operators
    Matrix operator*(const Matrix& A) const {return Matrix(M.cwiseProduct(A.M));}
    Matrix operator+(const Matrix& A) const {return Matrix(M + A.M);}
    Matrix operator-(const Matrix& A) const {return Matrix(M - A.M);}
    bool operator==(const Matrix& A) const {return M == A.M;}
    Matrix operator*(double a) const {return Matrix(a * M);}

    // other operators
    friend std::ostream& operator<<(std::ostream& os, const Matrix& A);
    double operator()(int i, int j) const {return M(i, j);}
    double& operator()(int i, int j) {return M(i, j);}

    // getters
    int rows() const {return M.rows();} int cols() const {return M.cols();}
    Block block(int i, int j, int m, int n) {return M.block(i, j, m, n);}
    EigenMat raw() const { return M; } EigenMat& raw() { return M; }
    Matrix left(int n) const {return Matrix(M.leftCols(n));}

    // mathematical operations
    std::tuple<Matrix, Matrix> eigensolve(const Matrix& A) const;
    Matrix dot(const Matrix& A) const {return Matrix(M * A.M);}
    Matrix inverse() const {return Matrix(M.inverse());}
    Matrix t() const {return Matrix(M.transpose());}
    double norm() const {return M.norm();}
    double sum() const {return M.sum();}

private:
    EigenMat M;
};

inline Matrix operator*(double a, const Matrix& A) {return A * a;}
