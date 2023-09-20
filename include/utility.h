#pragma once

#include "eigen.h"

namespace Utility {
    void SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn, const Vector& energy);
};
