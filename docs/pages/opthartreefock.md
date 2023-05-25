---
title: Hartree-Fock
parent: Usage
layout: default
nav_order: 2
---
{% include mathjax.html %}

# Hartree-Fock

Specify `"title" : "hf"` if you want to perform the Hartree-Fock calculation. You need to specify the molecule block before this calculation. Restricted and unrestricted versions are chosen mased on charge and multiplicity. Below are the available options.

## Main Block Options

{: .note-title }
> maxiter
>
> __Description__: Maximum number of iterations in the SCF cycle.<br>
> __Datatype__: int<br>
> __Default__: 100

{: .note-title }
> thresh
>
> __Description__: When difference in energy and RMS of density matrix in consecutive iterations reach this value the SCF is considered as converged.<br>
> __Datatype__: double<br>
> __Default__: 1e-8

{: .note-title }
> diis
>
> __Description__: Enables DIIS with default options.
>
> {: .note-title }
> > start
> >
> > __Description__: Iteration where the DIIS is turned on.<br>
> > __Datatype__: int<br>
> > __Default__: 3<br>
>
> {: .note-title }
> > keep
> >
> > __Description__: Number of previous Fock matrices to keep.<br>
> > __Datatype__: int<br>
> > __Default__: 5<br>
>
> {: .note-title }
> > damp
> >
> > __Description__: Damping parameter.<br>
> > __Datatype__: double<br>
> > __Default__: 0<br>

{: .note-title }
> print
>
> __Description__: Some printing options.
>
> {: .note-title }
> > iter
> >
> > __Description__: Whether or not to print the iterations.<br>
> > __Datatype__: bool<br>
> > __Default__: true<br>
>
> {: .note-title }
> > kinetic
> >
> > __Description__: Whether or not to print the one-electron kinetic matrix.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
>
> {: .note-title }
> > overlap
> >
> > __Description__: Whether or not to print the one-electron overlap matrix.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
>
> {: .note-title }
> > oneelec
> >
> > __Description__: Whether or not to print the one-electron nuclear attraction matrix.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
>
> {: .note-title }
> > orben
> >
> > __Description__: Whether or not to print the orbital energies.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
>
> {: .note-title }
> > mos
> >
> > __Description__: Whether or not to print the calculated matrix of coefficients.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
>
> {: .note-title }
> > density
> >
> > __Description__: Whether or not to print the density matrix.<br>
> > __Datatype__: bool<br>
> > __Default__: false<br>
