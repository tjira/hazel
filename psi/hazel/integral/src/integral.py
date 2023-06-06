import psi4

def integral(molecule):
    psi4.core.plugin("integral.so", psi4.core.Wavefunction.build(molecule, psi4.core.get_global_option("basis")))
