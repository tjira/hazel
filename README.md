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
Hazel is a collection of various electronic structure methods. It works similarly to the other available quantum programs. It is provided with an input file with system and calculation specification and prints the results to the terminal. Part of the software is also a molecular viewer.
</p>

## ✨ Features

Below are all the important features of Hazel divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock Method
* Møller–Plesset Perturbation Theory

### Convergence Accelerators

* Direct Inversion in the Iterative Subspace

### Additional Calculations

* Gradients for HF (analytical) and MP2 (numerical)
* Numerical hessians and frequency analysis for HF and MP2
* HF and MP2 steepest descent optimization

## 🛠️ Compilation

The software requires the [libint](https://github.com/evaleev/libint) library. You can compile it yourself following the instructions on the webpage (or better yet execute `./script/libint.sh` from the root directory). Before compilation make sure that you have Eigen library installed. You can do that on debian-based distros using the following command.

```bash
sudo apt install libeigen3-dev
```

Make sure that you add the libint headers and library to the environmental variables `CPLUS_INCLUDE_PATH` and `LIBRARY_PATH`. You can then configure the project by the following command.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release .
```

If the configuration finished without errors, compile the project by running the following command.

```bash
cmake --build build
```

After the compilation the bin folder will be created along with the executable.

## ⭐ Credits

* [libint](https://github.com/evaleev/libint) - High-performance library for computing Gaussian integrals in quantum mechanics.
* [libxc](https://gitlab.com/libxc/libxc) - Library of exchange-correlation functionals for density-functional theory.
* [eigen](https://gitlab.com/libeigen/eigen) - Eigen is a C++ template library for linear algebra.
* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
