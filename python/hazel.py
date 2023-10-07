import numpy as np
import scipy as sp

A2AU = 1.8897259886

SM2AN = {
    "H" : 1
}

class System:
    def __init__(self, path):
        self.content = "".join(open(path).readlines()[2:])[:-1]

    def atoms(self):
        return [SM2AN[symbol] for symbol in self.content.split()[0::4]]

    def coords(self):
        return np.array([self.content.split()[i::4] for i in range(1, 4)]).astype(float).T

    def repulsion(self):
        E = 0
        for i in range(len(self.atoms())):
            for j in range(i + 1, len(self.atoms())):
                distance = np.linalg.norm(A2AU * (self.coords()[j, :] - self.coords()[i, :]))
                E += self.atoms()[i] * self.atoms()[j] / distance
        return E

if __name__ == "__main__":
    system = System("molecule.xyz")

    T, V, S = np.loadtxt("T.mat"), np.loadtxt("V.mat"), np.loadtxt("S.mat")
    J, Hcore, nocc = np.loadtxt("J.mat").reshape(2 * T.shape), T + V, 1
    
    ERI, D, E = J - 0.5 * np.transpose(J, (0, 3, 2, 1)), np.zeros_like(S), 0
    A = sp.linalg.fractional_matrix_power(S, -0.5)

    for i in range(5):
        F = Hcore + np.einsum('pqrs,rs->pq', ERI, D, optimize=True)
        eps, C = np.linalg.eigh(A.dot(F).dot(A)); C = A.dot(C)

        D = 2 * C[:, :nocc] * C[:, :nocc].T
        E = 0.5 * (D * (Hcore + F)).sum()

    print(E + system.repulsion())
