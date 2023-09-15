#include "utility.h"

void Utility::SaveWavefunction(const std::string& fname, const CVector& r, const std::vector<std::vector<CVector>>& wfn) {
    int size = 0;
    for (size_t i = 0; i < wfn.size(); i++) {
        size = std::max(size, (int)wfn.at(i).size());
    }
    std::ofstream file(fname);
    for (int i = 0; i < size; i++) {
        file << std::fixed << std::setprecision(14) << "#        r1         ";
        for (int j = 0; j < wfn.size(); j++) {
            file << "     state" << (j < 10 ? "0" : "") << j << ".real         " << "state" << (j < 10 ? "0" : "") << j << ".imag    ";
        }
        for (int j = 0; j < r.size(); j++) {
            file << (!j ? "\n" : "") << std::setw(20) << r(j).real();
            for (int k = 0; k < wfn.size(); k++) {
                file << " " << std::setw(20) << wfn.at(k).at(std::min(i, (int)wfn.at(k).size() - 1))(j).real()
                << " " << std::setw(20) << wfn.at(k).at(std::min(i, (int)wfn.at(k).size() - 1))(j).imag();
            }
            file << "\n";
        }
    }
}