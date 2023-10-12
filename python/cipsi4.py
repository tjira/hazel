"""
A Psi4 input script to compute Full Configuration Interaction from a SCF reference

Requirements:
SciPy 0.13.0+, NumPy 1.7.2+

References:
Equations from [Szabo:1996]
"""

__authors__ = "Tianyuan Zhang"
__credits__ = ["Tianyuan Zhang", "Jeffrey B. Schriber", "Daniel G. A. Smith"]

__copyright__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__ = "BSD-3-Clause"
__date__ = "2017-05-26"

import time
import numpy as np
np.set_printoptions(precision=5, linewidth=200, suppress=True)
import psi4

np.set_printoptions(edgeitems=18)

# Check energy against psi4?
compare_psi4 = True

# Memory for Psi4 in GB
# psi4.core.set_memory(int(2e9), False)
psi4.core.set_output_file('output.dat', False)

# Memory for numpy in GB
numpy_memory = 2

mol = psi4.geometry("""
 H  0.02226951616055  0.00000000000000  0.00000000000000
 F  0.97773048383945  0.00000000000000  0.00000000000000
symmetry c1
""")


psi4.set_options({'basis': 'sto-3g',
                  'scf_type': 'pk',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

print('\nStarting SCF and integral build...')
t = time.time()

# First compute SCF energy using Psi4
scf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Grab data from wavfunction class
C = wfn.Ca()
print(np.asarray(C))
ndocc = wfn.doccpi()[0]
nmo = wfn.nmo()
# print(np.asarray(C))

# Compute size of Hamiltonian in GB
from scipy.special import comb
nDet = comb(nmo, ndocc)**2
H_Size = nDet**2 * 8e-9
print('\nSize of the Hamiltonian Matrix will be %4.2f GB.' % H_Size)
if H_Size > numpy_memory:
    clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                    limit of %4.2f GB." % (H_Size, numpy_memory))

# Integral generation from Psi4's MintsHelper
t = time.time()
mints = psi4.core.MintsHelper(wfn.basisset())
H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())

print('\nTotal time taken for ERI integrals: %.3f seconds.\n' % (time.time() - t))

#Make spin-orbital MO
print('Starting AO -> spin-orbital MO transformation...')
t = time.time()
MO = np.asarray(mints.mo_spin_eri(C, C))

# Update H, transform to MO basis and tile for alpha/beta spin
H = np.einsum('uj,vi,uv', C, C, H)
H = np.repeat(H, 2, axis=0)
H = np.repeat(H, 2, axis=1)

# Make H block diagonal
spin_ind = np.arange(H.shape[0], dtype=int) % 2
H *= (spin_ind.reshape(-1, 1) == spin_ind)
# print(MO[0, 1, :, :])

# ===============================================================================================================================================
# C_spin = np.repeat(C, 2, axis=0)
# C_spin = np.repeat(C_spin, 2, axis=1)
# J = np.asarray(mints.ao_eri())
# C_spin *= (spin_ind.reshape(-1, 1) == spin_ind)
# print(C_spin)
# C_spin = psi4.core.Matrix.from_array(C_spin)
# JMO = np.asarray(mints.mo_eri(C_spin,C_spin,C_spin,C_spin))

# print(ndocc)
# print(np.asarray(C))

# JMO = np.asarray(mints.mo_eri(C, C, C, C))
# JMO = np.repeat(JMO, 2, axis=0)
# JMO = np.repeat(JMO, 2, axis=1)
# JMO = np.repeat(JMO, 2, axis=2)
# JMO = np.repeat(JMO, 2, axis=3)
# spin_ind = np.arange(H.shape[0], dtype=np.int) % 2
# JMO *= (spin_ind.reshape(-1, 1, 1, 1) == spin_ind)
# print(np.asarray(wfn.Ca()))
# print(np.asarray(wfn.Cb()))
eps_a = np.asarray(wfn.epsilon_a())
eps_b = np.asarray(wfn.epsilon_b())
eps = np.append(eps_a, eps_b)

def spin_block_tei(I):
    identity = np.eye(2)
    I = np.kron(identity, I)
    return np.kron(identity, I.T)
 
J = spin_block_tei(np.asarray(mints.ao_eri()))
C_spin = np.block([
             [      wfn.Ca(),         np.zeros_like(C)],
             [np.zeros_like(C),          wfn.Cb()     ]])
C_spin = C_spin[:, eps.argsort()]
print(J.shape)

J = J.transpose(0, 2, 1, 3) - J.transpose(0, 2, 3, 1)
# JMO = np.einsum("pi,qa,pqrs,rj,sb->iajb", C_spin, C_spin, J, C_spin, C_spin, optimize=True)
# JMO = np.einsum('pQRS, pP -> PQRS',
#       np.einsum('pqRS, qQ -> pQRS',
#       np.einsum('pqrS, rR -> pqRS',
#       np.einsum('pqrs, sS -> pqrS', J, C_spin, optimize=True), C_spin, optimize=True), C_spin, optimize=True), C_spin, optimize=True)
JMO = np.einsum("pi,qa,pqrs,rj,sb->iajb", C_spin, C_spin, J, C_spin, C_spin, optimize=True)
# np.savetxt("JMO_spin.txt", JMO)



# print(JMO[0, 1, :, :])
# print("="*100)
# print(np.sum(JMO - MO))
# print((JMO == MO).all())
# ===============================================================================================================================================

print('..finished transformation in %.3f seconds.\n' % (time.time() - t))

from helper_CI import Determinant, HamiltonianGenerator
from itertools import combinations

print('Generating %d Full CI Determinants...' % (nDet))
t = time.time()
detList = []
for alpha in combinations(range(nmo), ndocc):
    for beta in combinations(range(nmo), ndocc):
        detList.append(Determinant(alphaObtList=alpha, betaObtList=beta))
# detList = detList[0].generateSingleExcitationsOfDet(nmo)

print('..finished generating determinants in %.3f seconds.\n' % (time.time() - t))

print('Generating Hamiltonian Matrix...')

t = time.time()
Hamiltonian_generator = HamiltonianGenerator(H, MO)
Hamiltonian_matrix = Hamiltonian_generator.generateMatrix(detList)
np.savetxt("psi.mat", Hamiltonian_matrix)


print('..finished generating Matrix in %.3f seconds.\n' % (time.time() - t))

print('Diagonalizing Hamiltonian Matrix...')

t = time.time()

e_fci, wavefunctions = np.linalg.eigh(Hamiltonian_matrix)
print('..finished diagonalization in %.3f seconds.\n' % (time.time() - t))
print(e_fci[0:10])

fci_mol_e = e_fci[0] + mol.nuclear_repulsion_energy()

print('# Determinants:     % 16d' % (len(detList)))

print('SCF energy:         % 16.10f' % (scf_e))
print('FCI correlation:    % 16.10f' % (fci_mol_e - scf_e))
print('Total FCI energy:   % 16.10f' % (fci_mol_e))

if compare_psi4:
    psi4.compare_values(psi4.energy('FCI'), fci_mol_e, 6, 'FCI Energy')














