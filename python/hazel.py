import numpy as np
import scipy as sp

def load():
    T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")
    J = np.loadtxt("J.mat").reshape(4 * [S.shape[1]]); return T, V, S, J

if __name__ == "__main__":
    # load the integrals and define the number of electrons
    [T, V, S, J], nel = load(), 2

    # define energy, number of occupied orbitals and convergence threshold
    E, Ep, n, thresh = 0, 1, nel // 2, 1e-6

    # define the core Hamiltonian, density and the coulomb and exchange matrices
    H, D = T + V, np.zeros_like(S)
    K = J.transpose(0, 3, 2, 1)

    while abs(E - Ep) > thresh:
        # build the Fock matrix
        F = H + np.einsum("ijkl,kl->ij", J - 0.5 * K, D)

        # solve the Fock equations
        eps, C = sp.linalg.eigh(F, S)

        # build the density and compute the energy
        D = 2 * np.einsum("ij,kj->ik", C[:, :n], C[:, :n])
        Ep, E = E, 0.5 * np.einsum("ij,ij->", D, H + F)

    # print the energy
    print(E)
