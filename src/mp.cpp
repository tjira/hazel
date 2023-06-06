#include "mp.h"

MP::MP(const System& system) : nocc(system.electrons / 2) {}

double MP::mp2(const Tensor<4>& Jmo, const Vector& eps) const {
    // define the energy
    double Ecorr = 0;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < Jmo.dimension(0); a++) {
                for (int b = nocc; b < Jmo.dimension(1); b++) {
                    Ecorr += Jmo(i, a, j, b) * (2 * Jmo(i, a, j, b) - Jmo(i, b, j, a)) / (eps(i) + eps(j) - eps(a) - eps(b));
                }
            }
        }
    }

    // return the energy
    return Ecorr;
}
