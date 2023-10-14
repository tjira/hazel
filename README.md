<h1 align="center">Quantum Hazel</h1>

<h4 align="center">
  <a href="https://github.com/tjira/hazel#%EF%B8%8F-compilation">Compilation</a>
  ·
  <a href="https://tjira.github.io/hazel/">Docs</a>
</h4>

<p align="center">
    <a href="https://github.com/tjira/hazel/pulse">
        <img src="https://img.shields.io/github/last-commit/tjira/hazel?logo=github&logoColor=white&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/hazel/blob/master/LICENSE.md">
        <img src="https://img.shields.io/github/license/tjira/hazel?logo=gitbook&logoColor=white&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/hazel/stargazers">
        <img src="https://img.shields.io/github/stars/tjira/hazel?logo=apachespark&logoColor=white&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/hazel">
        <img src="https://img.shields.io/github/languages/code-size/tjira/hazel?logo=databricks&logoColor=white&style=for-the-badge"/>
    </a>
    <br>
    <a href="https://github.com/tjira/hazel/releases/latest">
        <img src="https://img.shields.io/github/v/release/tjira/hazel?display_name=tag&logo=sharp&logoColor=white&style=for-the-badge"/>
    </a>
    <a href="https://github.com/tjira/hazel/releases/latest">
        <img src="https://img.shields.io/github/downloads/tjira/hazel/total?logo=markdown&logoColor=white&style=for-the-badge"/>
    </a>
</p>

<p align="center">
Hazel is a collection of various electronic structure methods. It works similarly to the other available quantum programs. It is provided with a geometry file, the calculation is specified int the command line and results are printed to the terminal.
</p>

## ✨ Features

Below are all the important features of Hazel divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock Method (RHF & UHF)
* Møller–Plesset Perturbation Theory
* Full Configuration Interaction

### Convergence Accelerators

* Direct Inversion in the Iterative Subspace

### Additional Calculations

* Numerical Hessians and Frequency Analysis for HF and MP2
* Gradients for HF (Analytical) and MP2 (Numerical)
* HF and RMP2 Steepest Descent Optimization
* Potential Energy Surface Scan
* Exact Quantum Dyamics
* Molecular Dynamics

## 🛠️ Compilation

### The Simple Way

The easiest cross-platform way to compile Hazel is using CMake with the following command.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DGOMP=ON -DSTANDALONE=ON .
```

If the configuration finished without errors, compile the project by running the following command.

```bash
cmake --build build
```

After the compilation the bin folder will be created along with the executable. Keep in mind that this way some of the features might be missing.

### The Proper Way

The proper way is only available on Linux. The software requires the [libint](https://github.com/evaleev/libint) library. Before the library compilation process, make sure you have [eigen](https://gitlab.com/libeigen/eigen) and [boost](https://github.com/boostorg/boost) installed. On debian-based distributions, you can do it with the following command.

```bash
sudo apt install libboost-all-dev libeigen3-dev
```

To compile the library execute `./script/libint.sh` from the project root directory. This command creates the *libint* folder with the compiled library. Now, we export the necessary environment variables.

```bash
export CPLUS_INCLUDE_PATH="$PWD/libint/install/include:$CPLUS_INCLUDE_PATH"
export LIBRARY_PATH="$PWD/libint/install/lib:$LIBRARY_PATH"
```

After this, the project configuration should finish without errors.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DGOMP=ON .
```

And we can build with the following command.

```bash
cmake --build build
```

After the compilation the bin folder will be created along with the executable.

## ⭐ Credits

* [libint](https://github.com/evaleev/libint) - High-performance library for computing Gaussian integrals in quantum mechanics.
* [libxc](https://gitlab.com/libxc/libxc) - Library of exchange-correlation functionals for density-functional theory.
* [eigen](https://gitlab.com/libeigen/eigen) - Eigen is a C++ template library for linear algebra.
* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
