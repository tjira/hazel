# Hazel

Hazel is a collection of various electronic structure methods. It works similarly to the other available quantum programs. It is provided an input file with system and calculation specification 
and prints the results to the terminal or redirected to a file. Part of the software is also a molecular viewer.

## Features

Below are all the important features of *Hazel* divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock

### Convergence Accelerators

* Direct Inversion in the Iterative Subspace

### Additional Calculations

* Numerical Gradients

## Compilation

To compile the software you will need to have cmake along with some basic libraries installed on your computer. To install all the necessary files on an Arch machine run the following command.

```bash
sudo pacman -S boost cmake eigen glm make nlohmann-json
```

Then you can prepare the makefile with the following command.

```bash
cmake -B build/release -DCMAKE_BUILD_TYPE=Release .
```

It will take a while since it's downloading some additional libraries. When the preparation is finished compile with the following command.

```bash
cmake --build build/release --parallel 4
```

The process can take a long time since the integral library *libint* is huge. After the compilation the *bin* folder will be created along with all the executables.

## Examples

To run a simple single point calculation, create a json file (for example *input.json*) with the following contents.

```json
{
    "system" : "molecule.xyz",
    "basis" : "STO-3G",
    "method" : {
        "name" : "HF",
        "diis" : {}
    }
}
```

You can then run the calculation with the `./hazel input.json` command.

## Acknowledgements

* [argparse](https://github.com/p-ranav/argparse) - Library for parsing command line arguments.
* [boost](https://github.com/boostorg/boost) - General library mainly used here for string formatting.
* [eigen](https://gitlab.com/libeigen/eigen) - Template library for linear algebra.
* [glad](https://gen.glad.sh/) - Multi-language OpenGL loader and generator.
* [glfw](https://github.com/glfw/glfw) - A multi-platform library for OpenGL window and input.
* [glm](https://github.com/g-truc/glm) - Mathematics for OpenGL.
* [imgui](https://github.com/ocornut/imgui) - Bloat-free GUI for C++ with minimal dependencies..
* [json](https://github.com/nlohmann/json) - JSON for modern C++.
* [libint](https://github.com/evaleev/libint) - High-performance library for computing Gaussian integrals in quantum mechanics.

