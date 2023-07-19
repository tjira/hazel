#include "ci.h"

Data CI::cid(bool) const {
    // define the output and number of occupied and virtual orbitals
    Data output = data; int nocc = data.system.electrons / 2; int nvirt = data.system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    output.ci.H = Matrix::Zero(1 * nocc * nvirt + 1, 1 * nocc * nvirt + 1);
    output.ci.H(0, 0) = data.hf.E - Integral::Repulsion(data.system);

    // fill the singlet double excitations
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
                    // fill the double-double elements
                    output.ci.H(i * nvirt + a + 1, j * nvirt + b + 1) = data.Jmo(nocc + b, nocc + a, nocc + a, nocc + b);

                    // fill the ground-double elements
                    output.ci.H(i * nvirt + a + 1, 0) = data.Jmo(nocc + a, i, i, nocc + a);
                    output.ci.H(0, j * nvirt + b + 1) = data.Jmo(nocc + b, j, j, nocc + b);
                }
            }
            // fill the diagonal elements
            output.ci.H(i * nvirt + a + 1, i * nvirt + a + 1) = output.ci.H(0, 0) + 2 * (data.hf.eps(nocc + a) - data.hf.eps(i)) + data.Jmo(i, i, i, i) + data.Jmo(nocc + a, nocc + a, nocc + a, nocc + a) - 4 * data.Jmo(nocc + a, nocc + a, i, i) + 2 * data.Jmo(nocc + a, i, i, nocc + a);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian, extract energies and return
    Eigen::SelfAdjointEigenSolver<Matrix> solver(output.ci.H); output.ci.C = solver.eigenvectors();
    output.ci.eig = solver.eigenvalues().array() + Integral::Repulsion(data.system);
    output.ci.Ecorr = output.ci.eig(0) - data.hf.E; return output;
}

Data CI::cis(bool) const {
    // define the output and number of occupied and virtual orbitals
    Data output = data; int nocc = data.system.electrons / 2; int nvirt = data.system.shells.nbf() - nocc;

    // create the CI Hamiltonian and fill the HF energy
    output.ci.H = Matrix::Zero(2 * nocc * nvirt + 1, 2 * nocc * nvirt + 1);
    output.ci.H(0, 0) = data.hf.E - Integral::Repulsion(data.system);

    // fill the singlet singles in CI Hamiltonian
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    output.ci.H(i * nvirt + a + 1, j * nvirt + b + 1) = 2 * data.Jmo(i, nocc + a, j, nocc + b) - data.Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            output.ci.H(i * nvirt + a + 1, i * nvirt + a + 1) += output.ci.H(0, 0) + data.hf.eps(nocc + a) - data.hf.eps(i);
        }
    }

    // fill the triplet singles in CI Hamiltonian
    for (size_t i = 0; i < nocc; i++) {
        for (size_t a = 0; a < nvirt; a++) {
            for (size_t j = 0; j < nocc; j++) {
                for (size_t b = 0; b < nvirt; b++) {
                    // fill the non-diagonal elements
                    output.ci.H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  j * nvirt + b + 1) -= data.Jmo(i, j, nocc + a, nocc + b);
                }
            }
            // fill the diagonal elements
            output.ci.H(nocc * nvirt + i * nvirt + a + 1,nocc * nvirt +  i * nvirt + a + 1) += output.ci.H(0, 0) + data.hf.eps(nocc + a) - data.hf.eps(i);
        }
    }

    // find the eigenvalues and eigenvectors of the CI Hamiltonian, extract energies and return
    Eigen::SelfAdjointEigenSolver<Matrix> solver(output.ci.H); output.ci.C = solver.eigenvectors();
    output.ci.eig = solver.eigenvalues().array() + Integral::Repulsion(data.system);
    output.ci.Ecorr = output.ci.eig(0) - data.hf.E; return output;
}
