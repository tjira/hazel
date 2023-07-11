#include "gradient.h"

template <class M>
Data Gradient<M>::get(bool print) const {
    // calculate the HF gradient
    if constexpr (std::is_same_v<HF, M>) {
        if (data.hf.grad.numerical) {
            auto efunc = [](Data data) {
                data.ints.S = Integral::Overlap(data.system);
                data.ints.T = Integral::Kinetic(data.system);
                data.ints.V = Integral::Nuclear(data.system);
                data.ints.J = Integral::Coulomb(data.system);
                return HF(data).scf(false);
            };
            return get(efunc, print);
        } else return getHF(print);

    // calculate the MP gradient
    } else if constexpr (std::is_same_v<MP, M>) {
        if (data.mp.grad.numerical) {
            auto efunc = [](Data data) {
                data.ints.S = Integral::Overlap(data.system), data.ints.T = Integral::Kinetic(data.system);
                data.ints.V = Integral::Nuclear(data.system), data.ints.J = Integral::Coulomb(data.system);
                data = HF(data).scf(false); data.intsmo.J = Transform::Coulomb(data.ints.J, data.hf.C);
                return MP(data).mp2();
            };
            return get(efunc, print);
        } else throw std::runtime_error("Analytical gradient for MP2 is not implemented");
    }
}

template <class M>
Data Gradient<M>::get(const std::function<Data(Data)>& efunc, bool print) const {
    // create the output data and define the gradient and step size
    Data output = data; Matrix G = Matrix::Zero(data.system.atoms.size(), 3); double step;

    // get the step value according to the method
    if constexpr (std::is_same_v<HF, M>) step = data.hf.grad.step;
    if constexpr (std::is_same_v<MP, M>) step = data.mp.grad.step;

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) shared(data, output) collapse(2)
    #endif
    for (int i = 0; i < G.rows(); i++) {
        for (int j = 0; j < G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices, temporary data.systems and integrals
            Matrix dirMinus(data.system.atoms.size(), 3); Data dataMinus = data;
            Matrix dirPlus(data.system.atoms.size(), 3); Data dataPlus = data;

            // fill the direction matrices
            dirMinus(i, j) -= step * A2BOHR; dirPlus(i, j) += step * A2BOHR;

            // move the systems
            dataMinus.system.move(dirMinus), dataPlus.system.move(dirPlus);
            dataMinus.hf.D = data.hf.D, dataPlus.hf.D = data.hf.D;

            // calculate the energies
            dataMinus = efunc(dataMinus); double energyMinus = dataMinus.hf.E;
            dataPlus = efunc(dataPlus); double energyPlus = dataPlus.hf.E;

            // add the correlation energy
            if constexpr (std::is_same_v<MP, M>) energyMinus += dataMinus.mp.Ecorr, energyPlus += dataPlus.mp.Ecorr;
                
            // calculate the derivative
            G(i, j) = BOHR2A * (energyPlus - energyMinus) / step / 2;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, G(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // assign the graient to the correct field
    if constexpr (std::is_same_v<HF, M>) output.hf.grad.G = G;
    if constexpr (std::is_same_v<MP, M>) output.mp.grad.G = G;

    // print empty line
    if (print) std::cout << std::endl;

    // return the gradient
    return output;
}

template <class M>
Data Gradient<M>::getHF(bool print) const {
    // extract the useful stuff from the calculated integrals and define all the contractio axes
    Tensor<3> dS1 = data.ints.dS.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dS.dimension(0), data.ints.dS.dimension(1), 3});
    Tensor<3> dT1 = data.ints.dT.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dT.dimension(0), data.ints.dT.dimension(1), 3});
    Tensor<3> dV1 = data.ints.dV.slice<Index<3>, Index<3>>({0, 0, 0}, {data.ints.dV.dimension(0), data.ints.dV.dimension(1), 3});
    Pair first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = data.system.electrons / 2;

    // define the density, weighed density, gradient matrix and ouptut
    Data output = data; output.hf.grad.G = Matrix::Zero(data.system.atoms.size(), 3);
    auto atom2shell = data.system.shells.atom2shell(data.system.atoms);
    Tensor<2> W(data.hf.C.rows(), data.hf.C.cols());

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * data.hf.C.leftCols(nocc).row(i).cwiseProduct(data.hf.C.leftCols(nocc).row(j)) * data.hf.eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI tensor
    Tensor<3> dERI = (data.ints.dJ - 0.5 * data.ints.dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(data.hf.D), Axes<2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < output.hf.grad.G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += data.system.shells.at(shell).size();

        // define the Hcore derivative and atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = data.ints.dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {data.hf.D.rows(), data.hf.D.cols(), 3});
        Eigen::array<Eigen::Index, 3> Soff = {si, 0, 0}, Sext = {ss, data.hf.D.cols(), 3};
        Eigen::array<Eigen::Index, 2> Doff = {si, 0}, Dext = {ss, data.hf.D.cols()};

        // fill the Hcore derivative
        dHcore.slice<Index<3>, Index<3>>({0, si, 0}, {data.hf.D.rows(), ss, 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext).shuffle(Index<3>{1, 0, 2});
        dHcore.slice<Index<3>, Index<3>>({si, 0, 0}, {ss, data.hf.D.cols(), 3}) += (dT1 + dV1).slice<Index<3>, Index<3>>(Soff, Sext);

        // contract the tensors and add them to the gradient
        output.hf.grad.G.row(i) += 2 * toVector(dERI.slice(Soff, Sext).contract(toTensor(data.hf.D).slice(Doff, Dext), Axes<2>{third, fourth}));
        output.hf.grad.G.row(i) -= 2 * toVector(dS1.slice(Soff, Sext).contract(W.slice(Doff, Dext), Axes<2>{third, fourth}));
        output.hf.grad.G.row(i) += toVector(dHcore.contract(toTensor(data.hf.D), Axes<2>{third, fourth}));
    }

    // add the nuclear repulsion contribution
    output.hf.grad.G += Integral::dRepulsion(data.system);

    // return the result
    return output;
}

template class Gradient<HF>;
template class Gradient<MP>;
