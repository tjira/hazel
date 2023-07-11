#include "hf.h"

Data HF::scf(bool print) const {
    // create the output data
    Data output = data;

    // create all the necessary matrices, calculate ERI and initialize DIIS
    Matrix H = data.ints.T + data.ints.V, F; int nocc = data.system.electrons / 2;
    Tensor<4> ERI = data.ints.J - 0.5 * data.ints.J.shuffle(Array<4>{0, 3, 2, 1});
    libint2::DIIS<Matrix> diis(data.hf.diis.start, data.hf.diis.keep);

    // specify the ERI contraction indices
    Eigen::IndexPair<int> first(2, 0), second(3, 1);

    // calculate the Fock matrix
    if (data.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(output.hf.D), Axes<2>{first, second}));
    else F = H + Integral::Coulomb(data.system, data.hf.D);

    // calculate the energy
    output.hf.E = 0.5 * output.hf.D.cwiseProduct(H + F).sum();

    // print the iteration header
    if (print) std::printf("\nITER       Eel [Eh]         |dE|     |dD|       TIME    \n");

    // perform the scf loop
    for (int i = 1; i <= data.hf.maxiter; i++) {
        // start the timer
        Timer::Timepoint start = Timer::Now();
        
        // calculate the Fock matrix
        if (data.ints.J.size()) F = H + toMatrix(ERI.contract(toTensor(output.hf.D), Axes<2>{first, second}));
        else F = H + Integral::Coulomb(data.system, output.hf.D);

        // exrapolate the fock matrix
        Matrix e = data.ints.S * output.hf.D * F - F * output.hf.D * data.ints.S;
        if (i > 1) diis.extrapolate(F, e);

        // solve the roothan equations
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver(F, data.ints.S);

        // exteract the eigenvalues and eigenvectors and save previous values of D and E
        output.hf.C = solver.eigenvectors(), output.hf.eps = solver.eigenvalues();
        Matrix Dp = output.hf.D; double Ep = output.hf.E;

        // calculate the new density and energy
        output.hf.D = 2 * output.hf.C.leftCols(nocc) * output.hf.C.leftCols(nocc).transpose();
        output.hf.E = 0.5 * output.hf.D.cwiseProduct(H + F).sum();

        // calculate the E and D errors and elapsed time
        double Eerr = std::abs(output.hf.E - Ep), Derr = (output.hf.D - Dp).norm();
        const char* elapsed = Timer::Format(Timer::Elapsed(start)).c_str();

        // print the iteration info line
        if (print) std::printf("%4d %20.14f %.2e %.2e %s %s\n", i, output.hf.E, Eerr, Derr, elapsed, i > data.hf.diis.start - 1 ? "DIIS" : "");

        // finish if covergence reached
        if (Eerr < data.hf.thresh && Derr < data.hf.thresh) break;
        else if (i == data.hf.maxiter) {
            throw std::runtime_error("Maximum number of iterations in SCF reached.");
        }
    }

    // add the nuclear repultion energy
    output.hf.E += Integral::Repulsion(data.system);

    // return the results
    return output;
}
