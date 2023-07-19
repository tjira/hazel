#include "ci.h"

CI::ResultsRestricted CI::cid(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the output and number of occupied and virtual orbitals
    int nocc = system.electrons / 2; int nvirt = system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    Matrix H(1 * nocc * nvirt + 1, 1 * nocc * nvirt + 1); H(0, 0) = ropt.rhfres.Eel;

    // fill the singlet double excitations
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
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

    // find the eigenvalues and eigenvectors of the CI Hamiltonian, extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eig = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eig(0) - ropt.rhfres.E;

    // return the results
    return {H, C, eig, Ecorr};
}

CI::ResultsRestricted CI::cis(const System& system, const Tensor<4>& Jmo, bool) const {
    // define the output and number of occupied and virtual orbitals
    int nocc = system.electrons / 2; int nvirt = system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    Matrix H(2 * nocc * nvirt + 1, 2 * nocc * nvirt + 1); H(0, 0) = ropt.rhfres.Eel;

    // fill the singlet singles in CI Hamiltonian
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    H(i * nvirt + a + 1, j * nvirt + b + 1) = 2 * Jmo(i, nocc + a, j, nocc + b) - Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            H(i * nvirt + a + 1, i * nvirt + a + 1) += H(0, 0) + ropt.rhfres.eps(nocc + a) - ropt.rhfres.eps(i);
        }
    }

    // fill the triplet singles in CI Hamiltonian
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  j * nvirt + b + 1) -= Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  i * nvirt + a + 1) += H(0, 0) + ropt.rhfres.eps(nocc + a) - ropt.rhfres.eps(i);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian, extract energies
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H); Matrix C = solver.eigenvectors();
    Vector eig = solver.eigenvalues().array() + ropt.rhfres.Enuc;
    double Ecorr = eig(0) - ropt.rhfres.E;

    // return the results
    return {H, C, eig, Ecorr};
}
