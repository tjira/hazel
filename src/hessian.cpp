#include "hessian.h"

template <class M>
std::tuple<Matrix, Vector> Hessian<M>::frequency(const System& system, bool) const {
    // create the hessian and frequencies
    Matrix H = get(system); Vector freq;

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

    // return the results
    return {H, freq};
}

template <class M>
Matrix Hessian<M>::get(const System& system, bool print) const {
    // calculate the HF hessian
    if constexpr (std::is_same_v<HF, M>) {
        if (ropt.numerical) {
            auto efunc = [this](System system) {
                system.ints.J = Integral::Coulomb(system), system.ints.S = Integral::Overlap(system);
                system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
                return HF(ropt.rhfopt).rscf(system, ropt.rhfres.D, false).E;
            };
            return get(system, efunc, print);
        } else throw std::runtime_error("ANALYTICAL HESSIAN FOR HF IS NOT IMPLEMENTED");

    // calculate the MP hessian
    } else if constexpr (std::is_same_v<MP, M>) {
        if (ropt.numerical) {
            auto efunc = [this](System system) {
                system.ints.J = Integral::Coulomb(system), system.ints.S = Integral::Overlap(system);
                system.ints.T = Integral::Kinetic(system), system.ints.V = Integral::Nuclear(system);
                HF::ResultsRestricted rhfres = HF(ropt.rhfopt).rscf(system, ropt.rhfres.D, false);
                Tensor<4> Jmo = Transform::Coulomb(system.ints.J, rhfres.C);
                return rhfres.E + MP({rhfres}).mp2(system, Jmo, false).Ecorr;
            };
            return get(system, efunc, print);
        } else throw std::runtime_error("ANALYTICAL HESSIAN FOR MP2 IS NOT IMPLEMENTED");
    }
}

template <class M>
Matrix Hessian<M>::get(const System& system, const std::function<double(System)>& efunc, bool print) const {
    // create the output data and define the gradient and step size
    Matrix H(3 * system.atoms.size(), 3 * system.atoms.size());

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
            Matrix dir1(system.atoms.size(), 3), dir2(system.atoms.size(), 3);

            // fill the direction matrices
            dir1(i / 3, i % 3) = ropt.step * A2BOHR; dir2(j / 3, j % 3) = ropt.step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);

            // calculate the energies
            double energyMinusMinus = efunc(sysMinusMinus), energyMinusPlus = efunc(sysMinusPlus);
            double energyPlusMinus = efunc(sysPlusMinus), energyPlusPlus = efunc(sysPlusPlus);

            // calculate the derivative
            H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / ropt.step / ropt.step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // print empty line
    if (print) std::cout << std::endl;

    // return the gradient
    return H;
}

template class Hessian<HF>;
template class Hessian<MP>;
