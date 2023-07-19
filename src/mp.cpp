#include "mp.h"

MP::ResultsRestricted MP::mp2(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the output, energy and nocc
    int nocc = system.electrons / 2; double Ecorr = 0;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < Jmo.dimension(0); a++) {
                for (int b = nocc; b < Jmo.dimension(1); b++) {
                    Ecorr += Jmo(i, a, j, b) * (2 * Jmo(i, a, j, b) - Jmo(i, b, j, a)) / (ropt.rhfres.eps(i) + ropt.rhfres.eps(j) - ropt.rhfres.eps(a) - ropt.rhfres.eps(b));
                }
            }
        }
    }

    // return the energy
    return {Ecorr};
}
