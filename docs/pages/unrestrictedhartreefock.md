---
title: Unrestricted Hartree–Fock
parent: Methods
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Unrestricted Hartree–Fock Method

Unrestricted Hartree–Fock (UHF) method is a method of approximation for the molecular wave function (or density matrix) and the energy of a quantum many-body system in a stationary state. In the unrestricted version electrons of different spins are treated separately, thus allowing of modelling molecules with different number of alpha and beta electrons. This page will assume that the reader is familiar with the RHF method.

## The Method

Since the one and two electrons integrals are not dependent on the density matrix they lok the same as in the restricted version. Let's start by defining $\alpha$ ad $\beta$ Fock matrix as

\begin{align}
F_{\mu\nu}^{\alpha}=H_{\mu\nu}^{core}+\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\alpha}(\mu\nu|\kappa\lambda)+\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\beta}(\mu\nu|\kappa\lambda)-\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\alpha}(\mu\lambda|\kappa\nu) \\\\\
F_{\mu\nu}^{\beta}=H_{\mu\nu}^{core}+\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\alpha}(\mu\nu|\kappa\lambda)+\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\beta}(\mu\nu|\kappa\lambda)-\frac{1}{2}\sum_{\kappa\lambda}D_{\kappa\lambda}^{\beta}(\mu\lambda|\kappa\nu).
\end{align}

The density matrices $\mathbf{D}^{\alpha}$ and $\mathbf{D}^{\beta}$ are formed the same way as in the restrited version with one small difference. One must take care when slicing the matrix of coefficients to form the density matrices. Number of columns in these matrices must match the number of corresponding electrons. With these density matrices, we solve two separate Fock equations, obtain expansion coefficients and orbital energies and reform the new density matrices. The total electron energy can be computed as

\begin{equation}
E=\frac{1}{4}\sum_{\mu\nu}D_{\mu\nu}^{\alpha}(H_{\mu\nu}^{core}+F_{\mu\nu}^{\alpha})+\frac{1}{4}\sum_{\mu\nu}D_{\mu\nu}^{\beta}(H_{\mu\nu}^{core}+F_{\mu\nu}^{\beta})+E_{nuc},
\end{equation}

where the nuclear energy is also the same as in the restricted version. 
