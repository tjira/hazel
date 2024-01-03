#pragma once

#include "eigen.h"

namespace Utility {
    inline std::string ToDblStr(double f) {std::stringstream ss; ss << std::fixed << std::setprecision(14) << std::setw(20) << f; return ss.str();}
    template <typename T> inline bool VectorContains(const std::vector<T>& v, const T& e) {return std::find(v.begin(), v.end(), e) != v.end();}
    void SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn, const Vector& energy);
    inline long Factorial(long n) {return n > 1 ? n * Factorial(n - 1) : 1;} std::vector<std::vector<int>> Combinations(int n, int k);
    inline std::string ToUpper(std::string s) {std::transform(s.begin(), s.end(), s.begin(), ::toupper); return s;}
    inline bool StringContains(const std::string& s, const std::string& e) {return s.find(e) != std::string::npos;}
    std::vector<int> Common(const std::vector<int>& first, const std::vector<int>& second);
    std::vector<int> Unique(const std::vector<int>& first, const std::vector<int>& second);
    Matrix Repeat(const Matrix& M, int count, int axis); 
};
