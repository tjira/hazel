#include "gradient.h"

template <class M>
Gradient<M>::Gradient(const Data& data) : data(data) {}

template <class M>
Data Gradient<M>::get(const System& system, bool print) const {
    // calculate the HF gradient
    if constexpr (std::is_same_v<HF, M>) {
        if (data.hf.grad.numerical) {
            auto efunc = [](System system, Data data) {
                system.ints.J = Integral::Coulomb(system);
                system.ints.S = Integral::Overlap(system);
                system.ints.T = Integral::Kinetic(system);
                system.ints.V = Integral::Nuclear(system);
                return HF(data).rscf(system, false);
            };
            return get(system, efunc, print);
        } else return getHF(system, print);

    // calculate the MP gradient
    } else if constexpr (std::is_same_v<MP, M>) {
        if (data.mp.grad.numerical) {
            auto efunc = [](System system, Data data) {
                system.ints.J = Integral::Coulomb(system);
                system.ints.S = Integral::Overlap(system);
                system.ints.T = Integral::Kinetic(system);
                system.ints.V = Integral::Nuclear(system);
                data = HF(data).rscf(system, false);
                HF::ResultsRestricted rhfres = {data.hf.C, data.hf.D, data.hf.eps, data.hf.E, data.hf.E - Integral::Repulsion(system), Integral::Repulsion(system)};
                Tensor<4> Jmo = Transform::Coulomb(system.ints.J, data.hf.C);
                data.mp.Ecorr = MP({rhfres}).mp2(system, Jmo, false).Ecorr;
                return data;
            };
            return get(system, efunc, print);
        } else throw std::runtime_error("ANALYTICAL GRADIENT FOR MP2 IS NOT IMPLEMENTED");
    }
}

template <class M>
Data Gradient<M>::get(const System& system, const std::function<Data(System, Data)>& efunc, bool print) const {
    // create the output data and define the gradient and step size
    Data output = data; Matrix G = Matrix::Zero(system.atoms.size(), 3); double step;

    // get the step value according to the method
    if constexpr (std::is_same_v<HF, M>) step = data.hf.grad.step;
    if constexpr (std::is_same_v<MP, M>) step = data.mp.grad.step;

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) collapse(2)
    #endif
    for (int i = 0; i < G.rows(); i++) {
        for (int j = 0; j < G.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices, temporary systems and integrals
            Matrix dirMinus(system.atoms.size(), 3); Data dataMinus = data; System sysMinus = system;
            Matrix dirPlus(system.atoms.size(), 3); Data dataPlus = data; System sysPlus = system;

            // fill the direction matrices
            dirMinus(i, j) -= step * A2BOHR; dirPlus(i, j) += step * A2BOHR;

            // move the systems
            dataMinus.hf.D = data.hf.D, dataPlus.hf.D = data.hf.D;
            sysMinus.move(dirMinus), sysPlus.move(dirPlus);

            // calculate the energies
            dataMinus = efunc(sysMinus, dataMinus); double energyMinus = dataMinus.hf.E;
            dataPlus = efunc(sysPlus, dataPlus); double energyPlus = dataPlus.hf.E;

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
Data Gradient<M>::getHF(const System& system, bool) const {
    // extract the useful stuff from the calculated integrals and define all the contractio axes
    Tensor<3> dS1 = system.dints.dS.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dS.dimension(0), system.dints.dS.dimension(1), 3});
    Tensor<3> dT1 = system.dints.dT.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dT.dimension(0), system.dints.dT.dimension(1), 3});
    Tensor<3> dV1 = system.dints.dV.slice<Index<3>, Index<3>>({0, 0, 0}, {system.dints.dV.dimension(0), system.dints.dV.dimension(1), 3});
    Pair first(2, 0), second(3, 1), third(0, 0), fourth(1, 1); int nocc = system.electrons / 2;

    // define the density, weighed density, gradient matrix and ouptut
    Data output = data; output.hf.grad.G = Matrix::Zero(system.atoms.size(), 3);
    auto atom2shell = system.shells.atom2shell(system.atoms);
    Tensor<2> W(data.hf.C.rows(), data.hf.C.cols());

    // fill the energy weighted density matrix
    for (int i = 0; i < W.dimension(0); i++) {
        for (int j = 0; j < W.dimension(1); j++) {
            W(i, j) = 2 * data.hf.C.leftCols(nocc).row(i).cwiseProduct(data.hf.C.leftCols(nocc).row(j)) * data.hf.eps.topRows(nocc);
        }
    }

    // calculate the derivative of the ERI tensor
    Tensor<3> dERI = (system.dints.dJ - 0.5 * system.dints.dJ.shuffle(Array<5>{0, 3, 2, 1, 4})).contract(toTensor(data.hf.D), Axes<2>{first, second});

    // for every gradient row (atom)
    for (int i = 0, si = 0, ss = 0; i < output.hf.grad.G.rows(); i++, si += ss, ss = 0) {
        // calculate number of shells for current atom
        for (long shell : atom2shell.at(i)) ss += system.shells.at(shell).size();

        // define the Hcore derivative and atomic slices for overlap tensor and density matrix
        Tensor<3> dHcore = system.dints.dV.slice<Index<3>, Index<3>>({0, 0, 6 + i * 3}, {data.hf.D.rows(), data.hf.D.cols(), 3});
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
    output.hf.grad.G += Integral::dRepulsion(system);

    // return the result
    return output;
}

template class Gradient<HF>;
template class Gradient<MP>;
