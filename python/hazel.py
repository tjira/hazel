import numpy as np
import itertools
# import psi4

import scipy as sp
import scipy.linalg
np.set_printoptions(edgeitems=18)
np.set_printoptions(linewidth=400)

BOHR2A, A2BOHR = 0.529177210903, 1 / 0.529177210903

SM2AN = {
    "H" : 1,
    "He" : 2,
    "C" : 6,
    "Ne" : 10,
    "O" : 8,
    "F" : 9,
    "Cl" : 17
}

# psi4.core.set_output_file('output.dat', False)
# mol = psi4.geometry("""
#  H  0.02226951616055  0.00000000000000  0.00000000000000
#  F  0.97773048383945  0.00000000000000  0.00000000000000
# symmetry c1
# """)
# psi4.set_options({'basis': 'sto-3g',
#                   'scf_type': 'pk',
#                   'e_convergence': 1e-8,
#                   'd_convergence': 1e-8})
# scf_e, wfn = psi4.energy('SCF', return_wfn=True)
# mints = psi4.core.MintsHelper(wfn.basisset())

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

def spin_block_tei(I):
    identity = np.eye(2)
    I = np.kron(identity, I)
    return np.kron(identity, I.T)

def align(arr1, arr2):
    arr1, arr2, swaps = list(arr1), list(arr2), 0
    for i in range(len(arr1)):
        if arr1[i] in arr2 and arr1[i] != arr2[i]:
            index = arr2.index(arr1[i])
            arr2[i], arr2[index] = arr2[index], arr2[i]
            swaps += 1
    return tuple(arr1), tuple(arr2), swaps

def aligndet(det1, det2):
    det1, det2 = list(det1), list(det2)
    det1[0], det2[0], swapsa = align(det1[0], det2[0])
    det1[1], det2[1], swapsb = align(det1[1], det2[1])
    return tuple(det1), tuple(det2), swapsa + swapsb

def occ(det, orb):
    counta, countb = 0, 0
    for a in det[0]:
        if a == orb: counta += 1
    for b in det[1]:
        if b == orb: countb += 1
    return counta, countb

def diff(det1, det2):
    if det1 == det2: return 0
    difference = 0
    for i in range(len(det1[0])):
        difference += int(det1[0][i] != det2[0][i]) + int(det1[1][i] != det2[1][i])
    # print(difflib.SequenceMatcher(None, det1[0], det2[0]).get_matching_blocks())
    return difference

def element(HMO, JMO, det1, det2):
    element = 0
    de1, det2, swaps = aligndet(det1, det2)
    solist1 = [2 * i for i in det1[0]] + [2 * i + 1 for i in det1[1]]
    solist2 = [2 * i for i in det2[0]] + [2 * i + 1 for i in det2[1]]
    common = list(set(solist1).intersection(solist2))
    # print(solist1)
    if diff(det1, det2) == 0:
        for so in solist1:
            element += HMO[so, so]
        for i in range(len(solist1)):
            for j in range(i + 1, len(solist1)):
                element += JMO[solist1[i], solist1[j], solist1[i], solist1[j]]
    if diff(det1, det2) == 1:
        exc = list()
        for i in range(len(det1[0])):
            if det1[0][i] > det2[0][i]: exc.append([2 * det2[0][i], 2 * det1[0][i]]) # a, r
            if det1[0][i] < det2[0][i]: exc.append([2 * det1[0][i], 2 * det2[0][i]]) # a, r
        for i in range(len(det1[1])):
            if det1[1][i] > det2[1][i]: exc.append([2 * det2[1][i] + 1, 2 * det1[1][i] + 1]) # a, r
            if det1[1][i] < det2[1][i]: exc.append([2 * det1[1][i] + 1, 2 * det2[1][i] + 1]) # a, r
        element += HMO[exc[0][0], exc[0][1]]
        for so in common:
            element += JMO[exc[0][0], so, exc[0][1], so]
        # element = 0
    if diff(det1, det2) == 2:
        exc = list()
        for i in range(len(det1[0])):
            if det1[0][i] > det2[0][i]: exc.append([2 * det2[0][i], 2 * det1[0][i]]) # a, r
            if det1[0][i] < det2[0][i]: exc.append([2 * det1[0][i], 2 * det2[0][i]]) # a, r
        for i in range(len(det1[1])):
            if det1[1][i] > det2[1][i]: exc.append([2 * det2[1][i] + 1, 2 * det1[1][i] + 1]) # a, r
            if det1[1][i] < det2[1][i]: exc.append([2 * det1[1][i] + 1, 2 * det2[1][i] + 1]) # a, r
        # print(exc[0][1], exc[1][1], exc[0][0], exc[1][0], (-1)**swaps * JMO[exc[0][1], exc[1][1], exc[0][0], exc[1][0]])
        element = JMO[exc[0][1], exc[1][1], exc[0][0], exc[1][0]]
        # element = 0
    return (-1)**swaps * element
    # return diff(det1, det2)

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
    print("FINAL HARTREE-FOCK ENERGY:", E + system.repulsion())

    # print(C)
    # align to psi4
    # C[:, 0] *= -1
    # C[:, -1] *= -1
    # C[:, [3, 4]] = C[:, [4, 3]]

    # physicists spin orbit notation for J
    J = spin_block_tei(J)
    C_spin = np.block([
                 [      C,         np.zeros_like(C)],
                 [np.zeros_like(C),          C     ]])
    C_spin = C_spin[:, np.append(eps, eps).argsort()]
    J = J.transpose(0, 2, 1, 3) - J.transpose(0, 2, 3, 1)

    # define the important matrices in MO basis
    HMO, JMO = Transform(C).oneelec(T + V), Transform(C_spin).twoelec(J)
    HMO = np.repeat(HMO, 2, axis=0)
    HMO = np.repeat(HMO, 2, axis=1)
    spin_ind = np.arange(HMO.shape[0], dtype=int) % 2
    HMO *= (spin_ind.reshape(-1, 1) == spin_ind)

    # generate all the singlet determinants 
    determinants = []
    for alpha in itertools.combinations(range(S.shape[0]), system.nocc()):
        for beta in itertools.combinations(range(S.shape[0]), system.nocc()):
            determinants.append((alpha, beta))
    # determinants2 = []
    # for i in range(1, len(determinants)):
    #     det1, det2, swaps = aligndet(determinants[0], determinants[i])
    #     if diff(det1, det2) == 1:
    #         determinants2.append(det2)
    # determinants = determinants2
    # print(determinants)
    # print(HMO)

    HCI = np.zeros(2 * [len(determinants)])

    # fill the matrix
    for i in range(HCI.shape[0]):
        for j in range(i + 1):
            HCI[i, j] = element(HMO, JMO, determinants[i], determinants[j])
            HCI[j, i] = HCI[i, j]

    # psihci = np.loadtxt("psi.mat")
    # print(np.sum(abs(psihci - HCI)))

    np.savetxt("hazel.mat", HCI)
    EPSCI, CCI = np.linalg.eigh(HCI)
    print(EPSCI[:10])
    print(EPSCI[0] + system.repulsion())
    
    # print(align((1, 2, 4, 3), (0, 1, 2, 3)))
    # print(JMO[0, 1, :, :])
    # print(ECI)
    # print(2 * HMO[0, 0] + 2 * H[1, 1] + 2 * H[2, 2] + 2 * H[3, 3] + 2 * H[4, 4] + )
    # print(list(itertools.permutations([*config])))



