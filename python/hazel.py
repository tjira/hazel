import numpy as np
import scipy as sp

BOHR2A, A2BOHR = 0.529177210903, 1 / 0.529177210903

SM2AN = {
    "H" : 1,
    "O" : 8
}

class Integral:
    def __init__(self, system):
        self.system = system

    def J(self):
        return np.loadtxt("J.mat").reshape(4 * tuple([int(np.sqrt(np.loadtxt("J.mat").shape[0]))]))

    def S(self):
        return np.loadtxt("S.mat")

    def T(self):
        return np.loadtxt("T.mat")

    def V(self):
        return np.loadtxt("V.mat")

class System:
    def __init__(self, path):
        self.content = "".join(open(path).readlines()[2:])[:-1]

    def atom(self, i):
        return np.array([SM2AN[symbol] for symbol in self.content.split()[0::4]])[i]

    def coords(self, i):
        return np.array([self.content.split()[i::4] for i in range(1, 4)]).astype(float).T[i, :]

    def distance(self, i, j):
        return np.linalg.norm(A2BOHR * (self.coords(j) - self.coords(i)))

    def nocc(self):
        return self.atom(range(self.size())).sum() // 2

    def repulsion(self):
        return sum([self.atom(i) * self.atom(j) / self.distance(i, j) for i in range(self.size()) for j in range(i + 1, self.size())])

    def size(self):
        return len(self.content.split()[0::4])



# the main function
if __name__ == "__main__":
    # load the system
    system = System("molecule.xyz")

    # calculate all the molecular integrals
    T, V, S, J = Integral(system).T(), Integral(system).V(), Integral(system).S(), Integral(system).J()

    # create the initial density guess and energy
    D, Dp, E, Ep = np.zeros_like(S), np.ones_like(S), 0, 1

    # invert the overlap matrix and define the electron repulsion tensor
    A, ERI = sp.linalg.fractional_matrix_power(S, -0.5), J - 0.5 * np.transpose(J, (0, 3, 2, 1))

    # start the scf loop
    while np.abs(E - Ep) > 1e-12 and np.linalg.norm(D - Dp) > 1e-12:
        # form the Fock matrix
        F = T + V + np.einsum("pqrs,rs->pq", ERI, D, optimize=True)

        # solve the Roothaan equations
        eps, C = np.linalg.eigh(A.dot(F).dot(A)); C = A.dot(C)

        # calculate the density matrix
        Dp = D; D = 2 * np.dot(C[:, :system.nocc()], C[:, :system.nocc()].T)

        # calculate the energy
        Ep = E; E = 0.5 * (D * (T + V + F)).sum()

    print(E + system.repulsion())
