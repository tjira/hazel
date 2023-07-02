#include "mp.h"

Data MP::mp2(bool) const {
    // define the output, energy and nocc
    Data output = data; output.mp.Ecorr = 0;
    int nocc = data.system.electrons / 2;

    // calculate the correlation energy
    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = nocc; a < data.intsmo.J.dimension(0); a++) {
                for (int b = nocc; b < data.intsmo.J.dimension(1); b++) {
                    output.mp.Ecorr += data.intsmo.J(i, a, j, b) * (2 * data.intsmo.J(i, a, j, b) - data.intsmo.J(i, b, j, a)) / (data.roothaan.eps(i) + data.roothaan.eps(j) - data.roothaan.eps(a) - data.roothaan.eps(b));
                }
            }
        }
    }

    // return the energy
    return output;
}

Data MP::Gradient::mp2(bool print) const {
    if (mp->data.mp.grad.numerical) return mp2Numerical(print);
    else return mp2Analytical(print);
}

Data MP::Gradient::mp2Analytical(bool) const {
    throw std::runtime_error("Analytical gradient for MP2 method is not implemented.");
}

Data MP::Gradient::mp2Numerical(bool print) const {
    // create the output data and define the gradient matrix
    Data output = mp->data; output.mp.grad.G = Matrix::Zero(mp->data.system.atoms.size(), 3);

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) shared(output)
    #endif
    for (int i = 0; i < output.mp.grad.G.rows(); i++) {
        for (int j = 0; j < output.mp.grad.G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices, temporary data.systems and integrals
            Matrix dirMinus(mp->data.system.atoms.size(), 3); Data dataMinus = mp->data;
            Matrix dirPlus(mp->data.system.atoms.size(), 3); Data dataPlus = mp->data;

            // fill the direction matrices
            dirMinus(i, j) -= mp->data.roothaan.grad.step * A2BOHR; dirPlus(i, j) += mp->data.roothaan.grad.step * A2BOHR;

            // move the systems
            dataMinus.roothaan.D = mp->data.roothaan.D, dataPlus.roothaan.D = mp->data.roothaan.D;
            dataMinus.system.move(dirMinus), dataPlus.system.move(dirPlus);

            // calculate all the integrals
            dataMinus.ints.S = Integral::Overlap(dataMinus.system), dataPlus.ints.S = Integral::Overlap(dataPlus.system);
            dataMinus.ints.T = Integral::Kinetic(dataMinus.system), dataPlus.ints.T = Integral::Kinetic(dataPlus.system);
            dataMinus.ints.V = Integral::Nuclear(dataMinus.system), dataPlus.ints.V = Integral::Nuclear(dataPlus.system);
            dataMinus.ints.J = Integral::Coulomb(dataMinus.system), dataPlus.ints.J = Integral::Coulomb(dataPlus.system);

            // calculate the HF energies
            dataMinus = Roothaan(dataMinus).scf(false);
            dataPlus = Roothaan(dataPlus).scf(false);

            // transform J to MO basis and calculate the correlation energy
            dataMinus.intsmo.J = Transform::Coulomb(dataMinus.ints.J, dataMinus.roothaan.C); dataMinus = MP(dataMinus).mp2(false);
            dataPlus.intsmo.J = Transform::Coulomb(dataPlus.ints.J, dataPlus.roothaan.C); dataPlus = MP(dataPlus).mp2(false);
                
            // calculate the derivative
            output.mp.grad.G(i, j) = BOHR2A * (dataPlus.roothaan.E + dataPlus.mp.Ecorr - dataMinus.roothaan.E - dataMinus.mp.Ecorr) / mp->data.roothaan.grad.step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, output.mp.grad.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // print empty line
    if (print) std::cout << std::endl;

    // return the gradient
    return output;
}

Data MP::Optimizer::mp2(bool print) const {
    // create the output and perform the SCF
    Data output = mp->data; output = Roothaan(output).scf(false);

    // transform the coulomb tensor
    output.intsmo.J = Transform::Coulomb(output.ints.J, output.roothaan.C);

    // calculate the MP2 energy and gradient
    output = MP(output).Gradient.mp2(false);
    output = MP(output).mp2(false);

    // print the header
    if (print) std::printf("ITER        E [Eh]         |GRAD|      TIME\n");

    // print the initial state info
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, output.roothaan.E + output.mp.Ecorr, output.mp.grad.G.norm(), "00:00:00.000");

    // move the data.system while gradient is big
    for (int i = 1; output.mp.grad.G.norm() > mp->data.mp.opt.thresh; i++) {
        // start the timer and move the system
        Timer::Timepoint start = Timer::Now();
        output.system.move(-output.mp.grad.G);

        // calculate integrals for HF method
        output.ints.S = Integral::Overlap(output.system);
        output.ints.T = Integral::Kinetic(output.system);
        output.ints.V = Integral::Nuclear(output.system);
        output.ints.J = Integral::Coulomb(output.system);

        // perform HF method
        output = Roothaan(output).scf(false);

        // transform the coulomb tensor
        output.intsmo.J = Transform::Coulomb(output.ints.J, output.roothaan.C);

        // calculate the MP2 energy and gradient
        output = MP(output).Gradient.mp2(false);
        output = MP(output).mp2(false);

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, output.roothaan.E + output.mp.Ecorr, output.mp.grad.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return results
    return output;
}
