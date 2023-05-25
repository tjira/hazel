---
title: Restricted Hartree–Fock
parent: Methods
layout: default
nav_order: 1
---
{% include mathjax.html %}

# Restricted Hartree–Fock Method

Restricted Hartree–Fock (RHF) method is a method of approximation for the molecular wave function (or density matrix) and the energy of a quantum many-body system in a stationary state. In he restricted version, electrons in the same orbitals are treated together.

## The Method

Let's start with defining the core Hamiltonian or in other words, the one-electron Hamiltonian. Core Hamiltonian is a part of the full Hamiltonian without the electron-electron repulsion. We can write it in index notation as

\begin{equation}
H_{\mu\nu}^{core}=T_{\mu\nu}+V_{\mu\nu},
\end{equation}

where $\mu$ and $\nu$ are indices of atomic orbitals, $T_{\mu\nu}$ is a kinetic energy matrix element and $V_{\mu\nu}$ is a potential energy matrix element. Both can be formally writen in a BraKet notation as

\begin{align}
T_{\mu\nu}&=\braket{\phi_{\mu}|\hat{T}|\phi_{\nu}} \\\\\
V_{\mu\nu}&=\braket{\phi_{\mu}|\hat{V}|\phi_{\nu}},
\end{align}

where $\phi_{\mu}$ is $\mu$-th atomic orbital. To practically calculate these matrix elements one first needs to introduce a basis set. When the basis set is introduced these integrals are handled by an integrated library libint. To continue we will also need to introduce the overlap matrix elements defined as

\begin{equation}
S_{\mu\nu}=\braket{\phi_{\mu}|\phi_{\nu}}
\end{equation}

and the two-electron coulomb repulsion tensor in the form

\begin{equation}
J_{\mu\nu\kappa\lambda}=\braket{\phi_{\mu}\phi_{\mu}|\hat{J}|\phi_{\kappa}\phi_{\lambda}}=(\mu\nu|\kappa\lambda).
\end{equation}

The notation $(\mu\nu\|\kappa\lambda)$ is used more frequently so I will use it here too. All these integrals can be calculated at the start of the calculation. The whole Hartree-Fock method revolves around solving the Roothan equations in the form

\begin{equation}
\mathbf{FC}=\mathbf{SC\varepsilon},
\end{equation}

where $\mathbf{F}$ is a Fock matrix, $\mathbf{C}$ is a matrix of orbital coefficients and $\mathbf{\varepsilon}$ will be orbital energies. The problem is that the Fock matrix in the form

\begin{equation}
F_{\mu\nu}=H_{\mu\nu}^{core}+\sum_{\kappa\lambda}D_{\kappa\lambda}\[(\mu\nu|\kappa\lambda)-\frac{1}{2}(\mu\lambda|\kappa\nu)\]
\end{equation}

depends on the density matrix $\mathbf{D}$, which we do not know. So we will have to solve it iteratively using a self consistent field (SCF) method. So, in the first iteration of the self consistent field we guess the Fock matrix and set it equal to the core Hamiltonian. Other first approximations can be used but this is probably the simplest. Then we solve the Roothan equations and calculate the density matrix as

\begin{equation}
D_{\mu\nu}=2\sum_{i}C_{\mu i}C_{\nu i}
\end{equation}

and the total energy of the system as

\begin{equation}
E=\frac{1}{2}\sum_{\mu\nu}D_{\mu\nu}(H_{\mu\nu}^{core}+F_{\mu\nu})+E_{nuc},
\end{equation}

where $E_{nuc}$ is nuclear repulsion energy which is dependent only on positions of nuclei so it can be precomputed before SCF as

\begin{equation}
E_{nuc}=\sum_{A}\sum_{B<A}\frac{Z_{A}Z_{B}}{R_{AB}}.
\end{equation}

We can now choose a convergence criterion as the difference in energies in consecutive iterations or maybe the Frobenius norm of the density matrix. If the difference in the criterion in two consecutive iterations is lower than some threshold (default $10^{-8}$ in hazel) we exit the loop otherwise we calculate the new Fock matrix with the new density matrix and repeat the process.

## The Gradient

If we perform the calculation as described above and get the density matrix $\mathbf{D}$ we can evaluate the nuclear energy gradient as

\begin{equation}
\frac{\partial E}{\partial X_A}=\sum_{\mu\nu}D_{\mu\nu}\frac{\partial H_{\mu\nu}^{core}}{\partial X_A}+2\sum_{\mu\nu\kappa\lambda}D_{\mu\nu}D_{\kappa\lambda}\frac{\partial(\mu\nu|\kappa\lambda)}{\partial X_A}-2\sum_{\mu\nu}W_{\mu\nu}\frac{\partial S_{\mu\nu}}{\partial X_A},
\end{equation}

where $\mathbf{W}$ is energy weighed density matrix defined as

\begin{equation}
W_{\mu\nu}=2\sum_{i}C_{\mu i}C_{\nu i}\varepsilon_i.
\end{equation}

All the integral derivatives are handled again by the libint library.
