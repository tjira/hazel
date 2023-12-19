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
Quantum Hazel, a dynamic collection of electronic structure methods, effortlessly transforms input geometry into quantum insights. Simply input a geometry file, specify calculations via the command line, and watch as results appear in your terminal.
</p>

## ✨ Features

Below are all the important features of Hazel divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock Method (RHF & UHF)
* Møller–Plesset Perturbation Theory
* Full Configuration Interaction

### Additional Calculations

* Gradients, Hessians and Frequency Analysis
* Steepest Descent Optimization
* Potential Energy Surface Scan
* Exact Quantum Dyamics
* Molecular Dynamics

### Convergence Accelerators

* Direct Inversion in the Iterative Subspace

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

The proper way is only available on Linux. The software requires the [libint](https://github.com/evaleev/libint) library. Before the library compilation process, make sure you have [eigen](https://gitlab.com/libeigen/eigen), [boost](https://github.com/boostorg/boost) an some basic X libraries installed. On debian-based distributions, you can do it with the following command.

```bash
sudo apt install libboost-all-dev libeigen3-dev libxcursor-dev libxi-dev libxinerama-dev libxrandr-dev
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

## 🖥️ Examples

This section show some basic calculation examples. You can view all the options of a specific command directly from the software by executing the calculation with added `-h` flag.

### Integral Calculation

This example calculates only the molecular integrals. To view or export them add the options from the following section to the `-p` (print) or `-e` (export) flag.

```bash
hazel -f molecule.xyz -b sto-3g ints
```

### Restricted Hartree-Fock Method

Basic HF calculation with no extra options. You can view all the options by adding the `-h` flag.

```bash
hazel -f molecule.xyz -b sto-3g rhf
```

### Post Hartree-Fock Methods

To execute some implemented post HF method add the corresponding keyword after the `rhf` keyword. To view available post HF methods add the `-h` flag after the `rhf` keyword. The example below shows the MP2 calculation.

```bash
hazel -f molecule.xyz -b sto-3g rhf mp2
```

### Quantum Gradients and Frequencies

Currently, the HF method is the only one with implemented analytical gradient. To calculate it execute the following command.

```bash
hazel -f molecule.xyz -b sto-3g rhf -g 0 1e-5
```

The 0 disables numerical gradient and 1e-5 is the step size for numerical gradient. The step is not ussed for analytical gradient, but you have to specify it if you plan to add another options after the `-g` flag. You can calculate numerical gradients with the post HF methods as follows.

```bash
hazel -f molecule.xyz -b sto-3g rhf fci -g 1 1e-5
```

The 1 enables numerical gradient and 1e-5 is the step size. The frequencies have to be calculated numerically. Below is an example to calculate it for the MP2 method.

```bash
hazel -f molecule.xyz -b sto-3g rhf mp2 -f 1 1e-5
```

### Molecular Optimization

To optimize a molecule using RHF method you can do the following. 

```bash
hazel -f molecule.xyz -b sto-3g opt -t 1e-8 rhf -g 1 1e-5
```

We also specified options for the gradient calculation. The 1e-8 is the gradient threshold for optimization.

### Molecular Dynamics

To perform simple MP2 molecular dynamics with default settings for the gradient execute the following.

```bash
hazel -f molecule.xyz -b sto-3g md rhf mp2
```

During the calculation the file trajectory.xyz will be created.

## 🖨️ Printing Options

In his section you find all the printing options implemented in Hazel. This also an overview of what data you, as a user, can extract from the software.

* `hazel`
    * `dist` - Distance matrix.
* `ints`
    * `dj` - Derivative of the coulomb integral in AO basis.
    * `ds` - Derivative of the overlap integral in AO basis.
    * `dt` - Derivaitve of the kinetic integral in AO basis.
    * `dv` - Derivative of the nuclear integral in AO basis.
    * `j` - Coulomb integral in AO basis (only if calculated).
    * `s` - Overlap integral in AO basis.
    * `t` - Kinetic integral in AO basis.
    * `v` - Nuclear integral in AO basis.
* `rhf`
    * `dist` - Optimized distance matrix (only for optimization).
    * `dj` - Derivative of the coulomb integral in AO basis (only if calculated).
    * `ds` - Derivative of the overlap integral in AO basis (only if calculated).
    * `dt` - Derivaitve of the kinetic integral in AO basis (only if calculated).
    * `dv` - Derivative of the nuclear integral in AO basis (only if calculated).
    * `c` - Matrix of coefficients in MO basis.
    * `d` - Density matrix in MO basis.
    * `eps` - Orbital energies.
    * `j` - Coulomb integral in AO basis (only if calculated).
    * `s` - Overlap integral in AO basis.
    * `t` - Kinetic integral in AO basis.
    * `v` - Nuclear integral in AO basis.
    * `pop` - Population analysis.
* `uhf`
    * `dist` - Optimized distance matrix (only for optimization).
    * `ca` - Matrix of coefficients for alpha electrons electrons in MO basis.
    * `cb` - Matrix of coefficients for beta electrons electrons in MO basis.
    * `da` - Density matrix for alpha electrons in MO basis.
    * `db` - Density matrix for beta electrons in MO basis.
    * `epsa` - Orbital energies for alpha electrons.
    * `epsb` - Orbital energies for beta electrons.
    * `j` - Coulomb integral in AO basis (only if calculated).
    * `s` - Overlap integral in AO basis.
    * `t` - Kinetic integral in AO basis.
    * `v` - Nuclear integral in AO basis.
* `rhf mp2`
    * `dist` - Optimized distance matrix (only for optimization).
    * `jmo` - Coulomb integral in MO basis.
* `rhf ci`, `rhf cis`, `rhf cid`, `rhf cisd`, `rhf fci`
    * `dist` - Optimized distance matrix (only for optimization).
    * `hms` - Hamiltonian in MS basis.
    * `jms` - Coulomb integral in MS basis.
    * `cih` - Configuration interaction Hamiltonian.
    * `cie` - Eigenvalues of CI Hamiltonian.
    * `cic` - Eigenvectors of CI Hamiltonian.

## ⭐ Credits

* [argparse](https://github.com/p-ranav/argparse) - Argument Parser for Modern C++.
* [glad](https://github.com/Dav1dde/glad) - Multi-Language Vulkan/GL/GLES/EGL/GLX/WGL Loader-Generator based on the official specs.
* [glfw](https://github.com/glfw/glfw) - A multi-platform library for OpenGL, OpenGL ES, Vulkan, window and input .
* [glm](https://github.com/g-truc/glm) - OpenGL Mathematics.
* [imgui](https://github.com/ocornut/imgui) - Bloat-free Graphical User Interface for C++ with minimal dependencies.
* [imguifiledialog](https://github.com/aiekick/ImGuiFileDialog) - File dialog for Dear ImGui.
* [json](https://github.com/nlohmann/json) - JSON for Modern C++.
* [stb](https://github.com/nothings/stb) - Single-file public domain libraries for C/C++.
