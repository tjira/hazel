#pragma once

#include "eigen.h"

namespace Utility {
    void SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn, const Vector& energy);
    inline bool VectorContains(const std::vector<std::string>& v, std::string e) {return std::find(v.begin(), v.end(), e) != v.end();}
    inline std::string ToUpper(std::string s) {std::transform(s.begin(), s.end(), s.begin(), ::toupper); return s;}
};
