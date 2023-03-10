# Hazel

Hazel is a collection of various electronic structure methods.
It consists of the *hazel* executable which is used to perform calculations and the *hview* executable which is a molecular viewer.

## Features

Below are all the important features of *Hazel* divided into categories.

### Quantum Mechanical Methods

* Hartree-Fock

### Convergence Accelerators

* Direct Inversion in the Iterative Subspace

### Additional Calculations

* Mulliken Analysis
* Numerical Gradients

## Compilation

Before compiling the software you will need some additional libraries. All the necessary libraries can be downloaded and set up with the `make libs` command.
The process can take up to an hour since it needs to download all the *boost* libraries and compile the *libint* library.
After the libraries are installed you can compile the software simply by running `make`. The *bin* folder will be created along with all the executables.

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

