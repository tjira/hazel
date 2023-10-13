#include "ci.h"

CI::ResultsRestricted CI::rcid(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print) const {
    // start the timer
    Timer::Timepoint start = Timer::Now();

    // generate all possible determinants
    std::vector<Determinant> dets = system.det().full();

    // define the remove function
    std::function<bool(Determinant)> remover = [dets](Determinant det) {
        int swaps; std::tie(det, swaps) = det.align(dets.at(0));
        int diff = dets.at(0).differences(det);
        return diff != 2;
    };

    // filter the determinants
    dets.erase(std::remove_if(dets.begin() + 1, dets.end(), remover), dets.end());

    // print the elapsed time
    std::cout << "\nGENERATED " << dets.size() << " DETERMINANTS: " << std::flush;
    std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // return the results
    return rsolve(dets, Hms, Jms, print);
}

CI::ResultsRestricted CI::rcis(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print) const {
    // start the timer
    Timer::Timepoint start = Timer::Now();

    // generate all possible determinants
    std::vector<Determinant> dets = system.det().full();

    // define the remove function
    std::function<bool(Determinant)> remover = [dets](Determinant det) {
        int swaps; std::tie(det, swaps) = det.align(dets.at(0));
        int diff = dets.at(0).differences(det);
        return diff != 1;
    };

    // filter the determinants
    dets.erase(std::remove_if(dets.begin() + 1, dets.end(), remover), dets.end());

    // print the elapsed time
    std::cout << "\nGENERATED " << dets.size() << " DETERMINANTS: " << std::flush;
    std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // return the results
    return rsolve(dets, Hms, Jms, print);
}

CI::ResultsRestricted CI::rcisd(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print) const {
    // start the timer
    Timer::Timepoint start = Timer::Now();

    // generate all possible determinants
    std::vector<Determinant> dets = system.det().full();

    // define the remove function
    std::function<bool(Determinant)> remover = [dets](Determinant det) {
        int swaps; std::tie(det, swaps) = det.align(dets.at(0));
        int diff = dets.at(0).differences(det);
        return diff != 1 && diff != 2;
    };

    // filter the determinants
    dets.erase(std::remove_if(dets.begin() + 1, dets.end(), remover), dets.end());

    // print the elapsed time
    std::cout << "\nGENERATED " << dets.size() << " DETERMINANTS: " << std::flush;
    std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // return the results
    return rsolve(dets, Hms, Jms, print);
}

CI::ResultsRestricted CI::rfci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool print) const {
    // start the timer
    Timer::Timepoint start = Timer::Now();

    // generate all possible determinants
    std::vector<Determinant> dets = system.det().full();

    // print the elapsed time
    std::cout << "\nGENERATED " << dets.size() << " DETERMINANTS: " << std::flush;
    std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // return the results
    return rsolve(dets, Hms, Jms, print);
}

CI::ResultsRestricted CI::rsolve(const std::vector<Determinant>& dets, const Matrix& Hms, const Tensor<4>& Jms, bool) const {
    // print the number of determinants
    std::cout << "\nFILLING RCI HAMILTONIAN: " << std::flush;

    // start the timer
    Timer::Timepoint start = Timer::Now();

    // fill the CI hamiltonian
    Matrix H(dets.size(), dets.size());
    for (int i = 0; i < H.rows(); i++) {
        for (int j = 0; j < i + 1; j++) {
            H(i, j) = dets.at(i).hamilton(dets.at(j), Hms, Jms); H(j, i) = H(i, j);
        }
    }

    // print the matrix creation time
    std::cout << Timer::Format(Timer::Elapsed(start)) << "\nFINDING THE EIGENVALUES: " << std::flush; start = Timer::Now();

    // find the eigenvalues and eigenvectors of the CI Hamiltonian and extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eps = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eps(0) - ropt.rhfres.E;

    // print the eigenproblem time
    std::cout << Timer::Format(Timer::Elapsed(start)) << std::endl;

    // return the results
    return {C, H, eps, Ecorr};
}
