#include "ci.h"

CI::ResultsRestricted CI::rcid(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the output and the number of occupied and virtual orbitals
    int nocc = system.electrons / 2; int nvirt = system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    Matrix H(1 * nocc * nvirt + 1, 1 * nocc * nvirt + 1); H(0, 0) = ropt.rhfres.Eel;

    // fill the singlet double excitations
    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nvirt; a++) {
            for (int j = 0; j < nocc; j++) {
                for (int b = 0; b < nvirt; b++) {
                    // fill the double-double elements
                    H(i * nvirt + a + 1, j * nvirt + b + 1) = Jmo(nocc + b, nocc + a, nocc + a, nocc + b);

                    // fill the ground-double elements
                    H(i * nvirt + a + 1, 0) = Jmo(nocc + a, i, i, nocc + a);
                    H(0, j * nvirt + b + 1) = Jmo(nocc + b, j, j, nocc + b);
                }
            }
            // fill the diagonal elements
            H(i * nvirt + a + 1, i * nvirt + a + 1) = H(0, 0) + 2 * (ropt.rhfres.eps(nocc + a) - ropt.rhfres.eps(i)) + Jmo(i, i, i, i) + Jmo(nocc + a, nocc + a, nocc + a, nocc + a) - 4 * Jmo(nocc + a, nocc + a, i, i) + 2 * Jmo(nocc + a, i, i, nocc + a);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian and extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eig = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eig(0) - ropt.rhfres.E;

    // return the results
    return {C, H, eig, Ecorr};
}

CI::ResultsRestricted CI::rcis(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the output and number of occupied and virtual orbitals
    int nocc = system.electrons / 2; int nvirt = system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    Matrix H(2 * nocc * nvirt + 1, 2 * nocc * nvirt + 1); H(0, 0) = ropt.rhfres.Eel;

    // fill the singlet singles in CI Hamiltonian
    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nvirt; a++) {
            for (int j = 0; j < nocc; j++) {
                for (int b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    H(i * nvirt + a + 1, j * nvirt + b + 1) = 2 * Jmo(i, nocc + a, j, nocc + b) - Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            H(i * nvirt + a + 1, i * nvirt + a + 1) += H(0, 0) + ropt.rhfres.eps(nocc + a) - ropt.rhfres.eps(i);
        }
    }

    // fill the triplet singles in CI Hamiltonian
    for (int i = 0; i < nocc; i++) {
        for (int a = 0; a < nvirt; a++) {
            for (int j = 0; j < nocc; j++) {
                for (int b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  j * nvirt + b + 1) -= Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  i * nvirt + a + 1) += H(0, 0) + ropt.rhfres.eps(nocc + a) - ropt.rhfres.eps(i);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian and extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eps = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eps(0) - ropt.rhfres.E;

    // return the results
    return {C, H, eps, Ecorr};
}

CI::ResultsRestricted CI::rfci(const System& system, const Matrix& Hms, const Tensor<4>& Jms, bool) const {
    // get all the possible determinants
    std::vector<Determinant> psis = system.det().full();

    // print the number of determinants
    std::cout << "\nCONSIDERING " << psis.size() << " DETERMINANTS" << std::endl;

    // antisymetrize the coulomb tensor and convert to physicist notation
    Tensor<4> Jmsap = Jms.shuffle(Array<4>{0, 2, 1, 3}) - Jms.shuffle(Array<4>{0, 2, 3, 1});

    // fill the CI hamiltonian
    Matrix H(psis.size(), psis.size());
    for (int i = 0; i < H.rows(); i++) {
        for (int j = 0; j < i + 1; j++) {
            H(i, j) = psis.at(i).hamilton(psis.at(j), Hms, Jmsap); H(j, i) = H(i, j);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian and extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eps = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eps(0) - ropt.rhfres.E;

    // return the results
    return {C, H, eps, Ecorr};
}
