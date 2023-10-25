#include "hessian.h"

Vector Hessian::frequency(const System& system, const Matrix& H, bool) const {
    // create the mass matrix
    Matrix MM(3 * system.atoms.size(), 3 * system.atoms.size());
    for (int i = 0; i < MM.rows(); i++) {
        MM(i, i) = std::sqrt(1 / masses.at(system.atoms.at(i / 3).atomic_number));
    }

    // calculate the frequencies in atomic units, exract them and return
    Eigen::EigenSolver<Matrix> solver(MM * H * MM); auto eval = solver.eigenvalues();
    Vector freq = (eval.cwiseSqrt() - eval.cwiseSqrt().imag()).real() * CFREQ;
    std::sort(freq.begin(), freq.end(), std::greater<>()); return freq;
}

Matrix Hessian::get(const System& system, const std::function<double(System)>& efunc, bool print) const {
    // initialize the hessian matrix
    Matrix H(3 * system.atoms.size(), 3 * system.atoms.size());

    // print the header
    if (print) std::printf("\n  ELEM      dE [Eh/Bohr]        TIME\n");

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
            dir1(i / 3, i % 3) = step * A2BOHR; dir2(j / 3, j % 3) = step * A2BOHR;

            // move the systems
            sysMinusMinus.move(-dir1 - dir2), sysMinusPlus.move(-dir1 + dir2);
            sysPlusMinus.move(dir1 - dir2), sysPlusPlus.move(dir1 + dir2);

            // calculate the energies
            double energyMinusMinus = efunc(sysMinusMinus), energyMinusPlus = efunc(sysMinusPlus);
            double energyPlusMinus = efunc(sysPlusMinus), energyPlusPlus = efunc(sysPlusPlus);

            // calculate and assign the derivative
            H(i, j) = BOHR2A * BOHR2A * (energyPlusPlus - energyMinusPlus - energyPlusMinus + energyMinusMinus) / step / step / 4;

            // print the iteration info
            if (print) std::printf("(%2d, %2d) %18.14f %s\n", i + 1, j + 1, H(i, j), Timer::Format(Timer::Elapsed(start)).c_str());
        }
    }

    // return the hessian
    return H;
}
