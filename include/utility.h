#pragma once

#include "eigen.h"

namespace Utility {
    void SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn, const Vector& energy);
    inline bool VectorContains(const std::vector<std::string>& v, std::string e) {return std::find(v.begin(), v.end(), e) != v.end();}
    inline int Factorial(int n) {return n > 1 ? n * Factorial(n - 1) : 1;} std::vector<std::vector<int>> Combinations(int n, int k);
    inline std::string ToUpper(std::string s) {std::transform(s.begin(), s.end(), s.begin(), ::toupper); return s;}
    Matrix Repeat(const Matrix& M, int count, int axis);
};
