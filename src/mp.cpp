#include "mp.h"

double MP::rmp2(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the correlation energy and nocc
    int nocc = system.electrons / 2; double Ecorr = 0;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < Jmo.dimension(0); a++) {
                for (int b = nocc; b < Jmo.dimension(1); b++) {
                    Ecorr += Jmo(i, a, j, b) * (2 * Jmo(i, a, j, b) - Jmo(i, b, j, a)) / (rhfres->eps(i) + rhfres->eps(j) - rhfres->eps(a) - rhfres->eps(b));
                }
            }
        }
    }

    // return the energy
    return Ecorr;
}
