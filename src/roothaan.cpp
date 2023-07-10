#include "roothaan.h"

Data Roothaan::frequency(bool print) const {
    // check if hessian is calculated
    if (!data.roothaan.freq.H.size()) throw std::runtime_error("You have not calculated the nuclear hessian matrix.");

    // create the output data
    Data output = data;

    // create the mass matrix
    Matrix M(3 * data.system.atoms.size(), 3 * data.system.atoms.size());
    for (int i = 0; i < M.rows(); i++) {
        M(i, i) = std::sqrt(1 / masses.at(data.system.atoms.at(i / 3).atomic_number));
    }

    // calculate the frequencies in atomic units
    Eigen::EigenSolver<Matrix> solver(M * output.roothaan.freq.H * M); auto eval = solver.eigenvalues();

    // extract the frequencies in cm-1
    output.roothaan.freq.freq = (eval.cwiseSqrt() - eval.cwiseSqrt().imag()).real() * CFREQ;

    // sort the frequencies
    std::sort(output.roothaan.freq.freq.begin(), output.roothaan.freq.freq.end(), std::greater<>());

    // return the results
    return output;
}

Data Roothaan::gradient(bool print) const {
    if (data.roothaan.grad.numerical) return gradientNumerical(print);
    else return gradientAnalytical(print);
}

Data Roothaan::gradientAnalytical(bool) const {
    // extract the useful stuff from the calculated integrals and define all the contractio axes
    Tensor<3> dS1 = data.ints.dS.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dS.dimension(0), data.ints.dS.dimension(1), 3});
    Tensor<3> dT1 = data.ints.dT.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dT.dimension(0), data.ints.dT.dimension(1), 3});
    Tensor<3> dV1 = data.ints.dV.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dV.dimension(0), data.ints.dV.dimension(1), 3});
    Pair first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = data.system.electrons / 2;

    // define the density, weighed density, gradient matrix and ouptut
    Data output = data; output.roothaan.grad.G = Matrix::Zero(data.system.atoms.size(), 3);
    auto atom2shell = data.system.shells.atom2shell(data.system.atoms);
    Tensor<2> W(data.roothaan.C.rows(), data.roothaan.C.cols());

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * data.roothaan.C.leftCols(nocc).row(i).cwiseProduct(data.roothaan.C.leftCols(nocc).row(j)) * data.roothaan.eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI tensor
    Tensor<3> dERI = (data.ints.dJ - 0.5 * data.ints.dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(data.roothaan.D), Axes<2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < output.roothaan.grad.G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += data.system.shells.at(shell).size();

        // define the Hcore derivative and atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = data.ints.dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {data.roothaan.D.rows(), data.roothaan.D.cols(), 3});
        Eigen::array<Eigen::Index, 3> Soff = {si, 0, 0}, Sext = {ss, data.roothaan.D.cols(), 3};
        Eigen::array<Eigen::Index, 2> Doff = {si, 0}, Dext = {ss, data.roothaan.D.cols()};

        // fill the Hcore derivative
        dHcore.slice<Index<3>, Index<3>>({0, si, 0}, {data.roothaan.D.rows(), ss, 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext).shuffle(Index<3>{1, 0, 2});
        dHcore.slice<Index<3>, Index<3>>({si, 0, 0}, {ss, data.roothaan.D.cols(), 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext);

        // contract the tensors and add them to the gradient
        output.roothaan.grad.G.row(i) += 2 * toVector(dERI.slice(Soff, Sext).contract(toTensor(data.roothaan.D).slice(Doff, Dext), Axes<2>{third, fourth}));
        output.roothaan.grad.G.row(i) -= 2 * toVector(dS1.slice(Soff, Sext).contract(W.slice(Doff, Dext), Axes<2>{third, fourth}));
        output.roothaan.grad.G.row(i) += toVector(dHcore.contract(toTensor(data.roothaan.D), Axes<2>{third, fourth}));
    }

    // add the nuclear repulsion contribution
    output.roothaan.grad.G += Integral::dRepulsion(data.system);

    // return the result
    return output;
}

Data Roothaan::gradientNumerical(bool print) const {
    // create the output data and define the gradient matrix
    Data output = data; output.roothaan.grad.G = Matrix::Zero(data.system.atoms.size(), 3);

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) shared(data, output) collapse(2)
    #endif
    for (int i = 0; i < output.roothaan.grad.G.rows(); i++) {
        for (int j = 0; j < output.roothaan.grad.G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices, temporary data.systems and integrals
            Matrix dirMinus(data.system.atoms.size(), 3); Data dataMinus = data;
            Matrix dirPlus(data.system.atoms.size(), 3); Data dataPlus = data;

            // fill the direction matrices
            dirMinus(i, j) -= data.roothaan.grad.step * A2BOHR; dirPlus(i, j) += data.roothaan.grad.step * A2BOHR;

            // move the systems
            dataMinus.roothaan.D = data.roothaan.D, dataPlus.roothaan.D = data.roothaan.D;
            dataMinus.system.move(dirMinus), dataPlus.system.move(dirPlus);

            // calculate all the integrals
            dataMinus.ints.S = Integral::Overlap(dataMinus.system), dataPlus.ints.S = Integral::Overlap(dataPlus.system);
            dataMinus.ints.T = Integral::Kinetic(dataMinus.system), dataPlus.ints.T = Integral::Kinetic(dataPlus.system);
            dataMinus.ints.V = Integral::Nuclear(dataMinus.system), dataPlus.ints.V = Integral::Nuclear(dataPlus.system);
            dataMinus.ints.J = Integral::Coulomb(dataMinus.system), dataPlus.ints.J = Integral::Coulomb(dataPlus.system);

            // calculate the energies
            dataMinus = Roothaan(dataMinus).scf(false);
            dataPlus = Roothaan(dataPlus).scf(false);
                
            // calculate the derivative
            output.roothaan.grad.G(i, j) = BOHR2A * (dataPlus.roothaan.E - dataMinus.roothaan.E) / data.roothaan.grad.step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, output.roothaan.grad.G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // print empty line
    if (print) std::cout << std::endl;

    // return the gradient
    return output;
}

Data Roothaan::hessian(bool print) const {
    if (data.roothaan.freq.numerical) return hessianNumerical(print);
    else return hessianAnalytical(print);
}

Data Roothaan::hessianAnalytical(bool) const {
    throw std::runtime_error("Analytical hessian for HF method is not implemented.");
}

Data Roothaan::hessianNumerical(bool print) const {
    // create the output data and define the gradient matrix
    Data output = data; output.roothaan.freq.H = Matrix::Zero(3 * data.system.atoms.size(), 3 * data.system.atoms.size());

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) shared(data, output) collapse(2)
    #endif
    for (int i = 0; i < output.roothaan.freq.H.rows(); i++) {
        for (int j = 0; j < output.roothaan.freq.H.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            Data dataMinusMinus = data, dataMinusPlus = data, dataPlusMinus = data, dataPlusPlus = data;
            Matrix dir1(data.system.atoms.size(), 3), dir2(data.system.atoms.size(), 3);

            // fill the direction matrices
            dir1(i / 3, i % 3) = data.roothaan.freq.step * A2BOHR; dir2(j / 3, j % 3) = data.roothaan.freq.step * A2BOHR;

            // move the systems
            dataMinusMinus.roothaan.D = data.roothaan.D, dataMinusPlus.roothaan.D = data.roothaan.D;
            dataPlusMinus.roothaan.D = data.roothaan.D, dataPlusPlus.roothaan.D = data.roothaan.D;
            dataMinusMinus.system.move(-dir1 - dir2), dataMinusPlus.system.move(-dir1 + dir2);
            dataPlusMinus.system.move(dir1 - dir2), dataPlusPlus.system.move(dir1 + dir2);

            // calculate all the integrals
            dataMinusMinus.ints.S = Integral::Overlap(dataMinusMinus.system), dataMinusPlus.ints.S = Integral::Overlap(dataMinusPlus.system);
            dataMinusMinus.ints.T = Integral::Kinetic(dataMinusMinus.system), dataMinusPlus.ints.T = Integral::Kinetic(dataMinusPlus.system);
            dataMinusMinus.ints.V = Integral::Nuclear(dataMinusMinus.system), dataMinusPlus.ints.V = Integral::Nuclear(dataMinusPlus.system);
            dataMinusMinus.ints.J = Integral::Coulomb(dataMinusMinus.system), dataMinusPlus.ints.J = Integral::Coulomb(dataMinusPlus.system);
            dataPlusMinus.ints.S = Integral::Overlap(dataPlusMinus.system), dataPlusPlus.ints.S = Integral::Overlap(dataPlusPlus.system);
            dataPlusMinus.ints.T = Integral::Kinetic(dataPlusMinus.system), dataPlusPlus.ints.T = Integral::Kinetic(dataPlusPlus.system);
            dataPlusMinus.ints.V = Integral::Nuclear(dataPlusMinus.system), dataPlusPlus.ints.V = Integral::Nuclear(dataPlusPlus.system);
            dataPlusMinus.ints.J = Integral::Coulomb(dataPlusMinus.system), dataPlusPlus.ints.J = Integral::Coulomb(dataPlusPlus.system);

            // calculate the energies
            dataMinusMinus = Roothaan(dataMinusMinus).scf(false);
            dataMinusPlus = Roothaan(dataMinusPlus).scf(false);
            dataPlusMinus = Roothaan(dataPlusMinus).scf(false);
            dataPlusPlus = Roothaan(dataPlusPlus).scf(false);

            // calculate the derivative
            output.roothaan.freq.H(i, j) = BOHR2A * BOHR2A * (dataPlusPlus.roothaan.E - dataMinusPlus.roothaan.E - dataPlusMinus.roothaan.E + dataMinusMinus.roothaan.E) / data.roothaan.freq.step / data.roothaan.freq.step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, output.roothaan.freq.H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // print empty line
    if (print) std::cout << std::endl;

    // return the output
    return output;

}

Data Roothaan::optimize(bool print) const {
    // create the output and perform the SCF and gradient calculation
    Data output = data; output = scf(false); output = Roothaan(output).gradient(false);

    // print the header
    if (print) std::printf("ITER        E [Eh]         |GRAD|      TIME\n");

    // print the initial state info
    if (print) std::printf("%4d %20.14f %.2e %s\n", 0, output.roothaan.E, output.roothaan.grad.G.norm(), "00:00:00.000");

    // move the data.system while gradient is big
    for (int i = 1; output.roothaan.grad.G.norm() > data.roothaan.opt.thresh; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();

        // move the system
        output.system.move(-output.roothaan.grad.G);

        // calculate integrals for HF method
        output.ints.S = Integral::Overlap(output.system);
        output.ints.T = Integral::Kinetic(output.system);
        output.ints.V = Integral::Nuclear(output.system);
        output.ints.J = Integral::Coulomb(output.system);

        // calculate integral derivatives for the gradient
        if (!data.roothaan.grad.numerical) {
            output.ints.dS = Integral::dOverlap(output.system);
            output.ints.dT = Integral::dKinetic(output.system);
            output.ints.dV = Integral::dNuclear(output.system);
            output.ints.dJ = Integral::dCoulomb(output.system);
        }

        // perform HF and calculate gradient
        output = Roothaan(output).scf(false); output = Roothaan(output).gradient(false);

        // print the iteration info
        if (print) std::printf("%4d %20.14f %.2e %s\n", i, output.roothaan.E, output.roothaan.grad.G.norm(), Timer::Format(Timer::Elapsed(start)).c_str());
    }

    // return results
    return output;
}

Data Roothaan::scf(bool print) const {
    // create the output data
    Data output = data;

    // create all the necessary matrices, calculate ERI and initialize DIIS
    Matrix H = data.ints.T + data.ints.V, F; int nocc = data.system.electrons / 2;
    Tensor<4> ERI = data.ints.J - 0.5 * data.ints.J.shuffle(Array<4>{0, 3, 2, 1});
    libint2::DIIS<Matrix> diis(data.roothaan.diis.start, data.roothaan.diis.keep);

    // specify the ERI contraction indices
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // calculate the Fock matrix
    if (data.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(output.roothaan.D), Axes<2>{first, second}));
    else F = H + Integral::Coulomb(data.system, data.roothaan.D);

    // calculate the energy
    output.roothaan.E = 0.5 * output.roothaan.D.cwiseProduct(H + F).sum();

    // print the iteration header
    if (print) std::printf("\nITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= data.roothaan.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();
        
        // calculate the Fock matrix
        if (data.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(output.roothaan.D), Axes<2>{first, second}));
        else F = H + Integral::Coulomb(data.system, output.roothaan.D);

        // exrapolate the fock matrix
        Matrix e = data.ints.S * output.roothaan.D * F - F * output.roothaan.D * data.ints.S;
        if (i > 1) diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, data.ints.S);

        // exteract the eigenvalues and eigenvectors and save previous values of D and E
        output.roothaan.C = solver.eigenvectors(), output.roothaan.eps = solver.eigenvalues();
        Matrix Dp = output.roothaan.D; double Ep = output.roothaan.E;

        // calculate the new density and energy
        output.roothaan.D = 2 * output.roothaan.C.leftCols(nocc) * output.roothaan.C.leftCols(nocc).transpose();
        output.roothaan.E = 0.5 * output.roothaan.D.cwiseProduct(H + F).sum();

        // calculate the E and D errors and elapsed time
        double Eerr = std::abs(output.roothaan.E - Ep), Derr = (output.roothaan.D - Dp).norm();
        const char* elapsed = Timer::Format(Timer::Elapsed(start)).c_str();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, output.roothaan.E, Eerr, Derr, elapsed, i > data.roothaan.diis.start - 1 ? "DIIS" : "");

        // finish if covergence reached
        if (Eerr < data.roothaan.thresh && Derr < data.roothaan.thresh) break;
        else if (i == data.roothaan.maxiter) {
            throw std::runtime_error("Maximum number of iterations in SCF reached.");
        }
    }

    // add the nuclear repultion energy
    output.roothaan.E += Integral::Repulsion(data.system);

    // return the results
    return output;
}
