#include "mp.h"

MP::MP(const Data& data) : data(data) {}

Data MP::mp2(bool) const {
    // define the output, energy and nocc
    Data output = data; output.mp.Ecorr = 0;
    int nocc = data.system.electrons / 2;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < data.Jmo.dimension(0); a++) {
                for (int b = nocc; b < data.Jmo.dimension(1); b++) {
                    output.mp.Ecorr += data.Jmo(i, a, j, b) * (2 * data.Jmo(i, a, j, b) - data.Jmo(i, b, j, a)) / (data.hf.eps(i) + data.hf.eps(j) - data.hf.eps(a) - data.hf.eps(b));
                }
            }
        }
    }

    // return the energy
    return output;
}
