import numpy as np
import scipy as sp

ATOM = {
    "H": 1,
    "O": 8,
}

def ints():
    T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")
    J = np.loadtxt("J.mat").reshape(4 * [S.shape[1]]); return T, V, S, J

def mol():
    with open("M.xyz") as M:
        nocc = sum([ATOM[A.split()[0]] for A in M.readlines()[2:]]) // 2
    return nocc

if __name__ == "__main__":
    # load the integrals and define the number of occupied orbitals
    [T, V, S, J], nocc = ints(), mol()

    # define energy and convergence threshold
    E, Ep, thresh = 0, 1, 1e-12

    # define the core Hamiltonian, density and the coulomb and exchange matrices
    H, D = T + V, np.zeros_like(S)
    K = J.transpose(0, 3, 2, 1)

    while abs(E - Ep) > thresh:
        # build the Fock matrix
        F = H + np.einsum("ijkl,ij", J - 0.5 * K, D)

        # solve the Fock equations
        eps, C = sp.linalg.eigh(F, S)

        # build the density from coefficients
        D = 2 * np.einsum("ij,kj", C[:, :nocc], C[:, :nocc])

        # calculate electron energy
        Ep, E = E, 0.5 * np.einsum("ij,ij", D, H + F)

    # print the energy
    print(E)
