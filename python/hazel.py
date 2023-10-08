import numpy as np
import itertools

import scipy as sp
import scipy.linalg

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

class Transform:
    def __init__(self, C):
        self.C = C

    def oneelec(self, M):
        return np.einsum("pi,pq,qj->ij", self.C, M, self.C, optimize=True)

    def twoelec(self, T):
        return np.einsum("pi,qa,pqrs,rj,sb->iajb", self.C, self.C, T, self.C, self.C, optimize=True)

def diff(det1, det2):
    if det1 == det2: return 0
    return 1000

def element(HMO, JMO, det1, det2):
    element = 0
    if diff(det1, det2) == 0:
        print(det1)
        for a in det1[0]:
            element += HMO[a, a]
        for b in det1[1]:
            element += HMO[b, b]
        for i in range(len(det1[0])):
            for j in range(i, len(det1[0])):
                element += JMO[det1[0][i], det1[0][i], det1[0][j], det1[0][j]] - JMO[det1[0][i], det1[0][j], det1[0][j], det1[0][i]]
        for i in range(len(det1[1])):
            for j in range(i, len(det1[1])):
                element += JMO[det1[1][i], det1[1][i], det1[1][j], det1[1][j]] - JMO[det1[1][i], det1[1][j], det1[1][j], det1[1][i]]
        for i in range(len(det1[0])):
            for j in range(i, len(det1[1])):
                element += JMO[det1[0][i], det1[0][i], det1[1][j], det1[1][j]]



        # for i in range(len(det1)):
        #     for j in range(i + 1, len(det2)):


        # print(list(itertools.product(det1, det2)))
        return element
    return element

# the main function
if __name__ == "__main__":
    # load the system
    system = System("molecule.xyz")

    # calculate all the molecular integrals
    T, V, S, J = Integral(system).T(), Integral(system).V(), Integral(system).S(), Integral(system).J()

    # create the initial density guess and energy
    D, Dp, E, Ep = np.zeros_like(S), np.ones_like(S), 0, 1

    # invert the overlap matrix and define the exchange tensor
    A, K = scipy.linalg.fractional_matrix_power(S, -0.5), np.transpose(J, (0, 3, 2, 1))

    # start the scf loop
    while np.abs(E - Ep) > 1e-12 and np.linalg.norm(D - Dp) > 1e-12:
        # form the Fock matrix
        F = T + V + np.einsum("pqrs,rs->pq", J - 0.5 * K, D, optimize=True)

        # solve the Roothaan equations
        eps, C = np.linalg.eigh(A.dot(F).dot(A)); C = A.dot(C)

        # calculate the density matrix
        Dp = D; D = 2 * np.dot(C[:, :system.nocc()], C[:, :system.nocc()].T)

        # calculate the energy
        Ep = E; E = 0.5 * (D * (T + V + F)).sum()

    # print the final energy
    print("FINAL HARTREE-FOCK ENERGY:", E)

    # define the important matrices in MO basis
    HMO, JMO, KMO = Transform(C).oneelec(T + V), Transform(C).twoelec(J), Transform(C).twoelec(K)
    
    # generate all the singlet determinants 
    determinants = []
    for alpha in itertools.combinations(range(S.shape[0]), system.nocc()):
        for beta in itertools.combinations(range(S.shape[0]), system.nocc()):
            determinants.append((alpha, beta))

    HCI = np.zeros(2 * [len(determinants)])

    # fill the matrix
    for i in range(HCI.shape[0]):
        for j in range(i + 1):
            HCI[i, j] = element(HMO, JMO, determinants[i], determinants[j])
            HCI[j, i] = element(HMO, JMO, determinants[i], determinants[j])

    print(HCI)
    # print(JMO[0, 1, :, :])
    print(HMO[0, 0] + HMO[1, 1] + JMO[0, 0, 1, 1])
    # print(2 * HMO[0, 0] + 2 * H[1, 1] + 2 * H[2, 2] + 2 * H[3, 3] + 2 * H[4, 4] + )
    # print(list(itertools.permutations([*config])))



