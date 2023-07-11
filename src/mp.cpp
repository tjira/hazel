#include "mp.h"

Data MP::mp2(bool) const {
    // check if the coulomb tensor was transformed
    if (!data.intsmo.J.size()) throw std::runtime_error("You have not transformed the coulomb tensor to the MO basis.");

    // define the output, energy and nocc
    Data output = data; output.mp.Ecorr = 0;
    int nocc = data.system.electrons / 2;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < data.intsmo.J.dimension(0); a++) {
                for (int b = nocc; b < data.intsmo.J.dimension(1); b++) {
                    output.mp.Ecorr += data.intsmo.J(i, a, j, b) * (2 * data.intsmo.J(i, a, j, b) - data.intsmo.J(i, b, j, a)) / (data.hf.eps(i) + data.hf.eps(j) - data.hf.eps(a) - data.hf.eps(b));
                }
            }
        }
    }

    // return the energy
    return output;
}
