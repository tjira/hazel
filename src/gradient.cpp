#include "gradient.h"

Matrix Gradient::get(const System& system, const std::function<double(System)>& efunc, bool print) const {
    // define the gradient matrix
    Matrix G(system.atoms.size(), 3);

    // print the header
    if (print) std::printf("\n  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) collapse(2)
    #endif
    for (int i = 0; i < G.rows(); i++) {
        for (int j = 0; j < G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            Matrix dirMinus(system.atoms.size(), 3); System sysMinus = system;
            Matrix dirPlus(system.atoms.size(), 3); System sysPlus = system;

            // fill the direction matrices
            dirMinus(i, j) -= opt.step * A2BOHR; dirPlus(i, j) += opt.step * A2BOHR;

            // move the systems
            sysMinus.move(dirMinus), sysPlus.move(dirPlus);
                
            // calculate and assign the derivative
            G(i, j) = BOHR2A * (efunc(sysPlus) - efunc(sysMinus)) / opt.step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the gradient
    return G;
}

Matrix Gradient::get(const System& system, const HF::ResultsRestricted& rhfres, bool) const {
    // extract the useful stuff from the calculated integrals and define all the contraction axes
    Tensor<3> dS1 = system.dints.dS.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dS.dimension(0), system.dints.dS.dimension(1), 3});
    Tensor<3> dT1 = system.dints.dT.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dT.dimension(0), system.dints.dT.dimension(1), 3});
    Tensor<3> dV1 = system.dints.dV.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dV.dimension(0), system.dints.dV.dimension(1), 3});
    Eigen::IndexPair<int> first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = system.electrons / 2;

    // define the weighted density, gradient matrix and atomic shell map
    auto atom2shell = system.shells.atom2shell(system.atoms);
    Tensor<2> W(rhfres.C.rows(), rhfres.C.cols());
    Matrix G(system.atoms.size(), 3);

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * rhfres.C.leftCols(nocc).row(i).cwiseProduct(rhfres.C.leftCols(nocc).row(j)) * rhfres.eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI tensor
    Tensor<3> dERI = (system.dints.dJ - 0.5 * system.dints.dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(rhfres.D), Axes<2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.shells.at(shell).size();

        // define the core Hamiltonian derivative, atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = system.dints.dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {rhfres.D.rows(), rhfres.D.cols(), 3});
        Index<3> Soff = {si, 0, 0}, Sext = {ss, rhfres.D.cols(), 3}; Index<2> Doff = {si, 0}, Dext = {ss, rhfres.D.cols()};

        // fill the core Hamiltonian derivative
        dHcore.slice<Index<3>, Index<3>>({0, si, 0}, {rhfres.D.rows(), ss, 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext).shuffle(Index<3>{1, 0, 2});
        dHcore.slice<Index<3>, Index<3>>({si, 0, 0}, {ss, rhfres.D.cols(), 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext);

        // contract the tensors and add them to the gradient
        G.row(i) += 2 * toVector(dERI.slice(Soff, Sext).contract(toTensor(rhfres.D).slice(Doff, Dext), Axes<2>{third, fourth}));
        G.row(i) -= 2 * toVector(dS1.slice(Soff, Sext).contract(W.slice(Doff, Dext), Axes<2>{third, fourth}));
        G.row(i) += toVector(dHcore.contract(toTensor(rhfres.D), Axes<2>{third, fourth}));
    }

    // add the nuclear contribution and return
    return G;
}
