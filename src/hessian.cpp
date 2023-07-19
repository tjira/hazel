#include "hessian.h"

template <class M>
Data Hessian<M>::frequency(const System& system, bool) const {
    // declare the output, hessian and frequencies
    Data output = data; Matrix H; Vector freq;

    // assign the hessian according to the method
    if constexpr (std::is_same_v<HF, M>) H = data.hf.freq.H;
    if constexpr (std::is_same_v<MP, M>) H = data.mp.freq.H;

    // check if hessian is calculated
    if (!H.size()) throw std::runtime_error("YOU HAVE NOT CALCULATED THE NUCLEAR HESSIAN MATRIX");

    // create the mass matrix
    Matrix MM(3 * system.atoms.size(), 3 * system.atoms.size());
    for (int i = 0; i < MM.rows(); i++) {
        MM(i, i) = std::sqrt(1 / masses.at(system.atoms.at(i / 3).atomic_number));
    }

    // calculate the frequencies in atomic units and exract them
    Eigen::EigenSolver<Matrix> solver(MM * H * MM); auto eval = solver.eigenvalues();
    freq = (eval.cwiseSqrt() - eval.cwiseSqrt().imag()).real() * CFREQ;
    std::sort(freq.begin(), freq.end(), std::greater<>());

    // assign the graient to the correct field
    if constexpr (std::is_same_v<HF, M>) output.hf.freq.freq = freq;
    if constexpr (std::is_same_v<MP, M>) output.mp.freq.freq = freq;

    // return the results
    return output;
}

template <class M>
Data Hessian<M>::get(const System& system, bool print) const {
    // calculate the HF hessian
    if constexpr (std::is_same_v<HF, M>) {
        if (data.hf.freq.numerical) {
            auto efunc = [](System system, Data data) {
                system.ints.J = Integral::Coulomb(system);
                system.ints.S = Integral::Overlap(system);
                system.ints.T = Integral::Kinetic(system);
                system.ints.V = Integral::Nuclear(system);
                return HF(data).rscf(system, false);
            };
            return get(system, efunc, print);
        } else throw std::runtime_error("ANALYTICAL HESSIAN FOR HF IS NOT IMPLEMENTED");

    // calculate the MP hessian
    } else if constexpr (std::is_same_v<MP, M>) {
        if (data.mp.freq.numerical) {
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
        } else throw std::runtime_error("ANALYTICAL HESSIAN FOR MP2 IS NOT IMPLEMENTED");
    }
}

template <class M>
Data Hessian<M>::get(const System& system, const std::function<Data(System, Data)>& efunc, bool print) const {
    // create the output data and define the gradient and step size
    Data output = data; Matrix H(3 * system.atoms.size(), 3 * system.atoms.size()); double step;

    // get the step value according to the method
    if constexpr (std::is_same_v<HF, M>) step = data.hf.freq.step;
    if constexpr (std::is_same_v<MP, M>) step = data.mp.freq.step;

    // print the header
    if (print) std::printf("  ELEM      dE [Eh/Bohr]        TIME\n");

    // fill the gradient
    #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthread) collapse(2)
    #endif
    for (int i = 0; i < H.rows(); i++) {
        for (int j = 0; j < H.cols(); j++) {
            // start the timer
            Timer::Timepoint start = Timer::Now();

            // define the direction matrices and temporary systems
            System sysMinusMinus = system, sysMinusPlus = system, sysPlusMinus = system, sysPlusPlus = system;
            Data dataMinusMinus = data, dataMinusPlus = data, dataPlusMinus = data, dataPlusPlus = data;
            Matrix dir1(system.atoms.size(), 3), dir2(system.atoms.size(), 3);

            // fill the direction matrices
            dir1(i / 3, i % 3) = step * A2BOHR; dir2(j / 3, j % 3) = step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);
            dataMinusMinus.hf.D = data.hf.D, dataMinusPlus.hf.D = data.hf.D;
            dataPlusMinus.hf.D = data.hf.D, dataPlusPlus.hf.D = data.hf.D;

            // calculate the energies
            dataMinusMinus = efunc(sysMinusMinus, dataMinusMinus); double energyMinusMinus = dataMinusMinus.hf.E;
            dataMinusPlus = efunc(sysMinusPlus, dataMinusPlus); double energyMinusPlus = dataMinusPlus.hf.E;
            dataPlusMinus = efunc(sysPlusMinus, dataPlusMinus); double energyPlusMinus = dataPlusMinus.hf.E;
            dataPlusPlus = efunc(sysPlusPlus, dataPlusPlus); double energyPlusPlus = dataPlusPlus.hf.E;

            // add the correlation energy if needed
            if constexpr (std::is_same_v<MP, M>) {
                energyMinusMinus += dataMinusMinus.mp.Ecorr, energyMinusPlus += dataMinusPlus.mp.Ecorr;
                energyPlusMinus += dataPlusMinus.mp.Ecorr, energyPlusPlus += dataPlusPlus.mp.Ecorr;
            }

            // calculate the derivative
            H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / step / step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // assign the graient to the correct field
    if constexpr (std::is_same_v<HF, M>) output.hf.freq.H = H;
    if constexpr (std::is_same_v<MP, M>) output.mp.freq.H = H;

    // print empty line
    if (print) std::cout << std::endl;

    // return the gradient
    return output;
}

template class Hessian<HF>;
template class Hessian<MP>;
