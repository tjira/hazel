#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"

namespace psi{
    namespace integral {
        extern "C" PSI_API SharedWavefunction integral(SharedWavefunction wfn, Options& options);
        extern "C" PSI_API int read_options(std::string name, Options &options);
    }
}

int psi::integral::read_options(std::string name, Options &options) {
    if (name == "INTEGRAL"|| options.read_globals()) {
        options.add_bool("NUCLEAR", false);
        options.add_bool("OVERLAP", false);
        options.add_bool("KINETIC", false);
        options.add_bool("ERI", false);
    }
    return true;
}

psi::SharedWavefunction psi::integral::integral(SharedWavefunction wfn, Options& options) {
    // define the integral engine
    MintsHelper mints(MintsHelper(wfn->basisset(), options, 0));

    // print the integrals
    if (options.get_bool("NUCLEAR")) mints.ao_potential()->print();
    if (options.get_bool("OVERLAP")) mints.ao_overlap()->print();
    if (options.get_bool("KINETIC")) mints.ao_kinetic()->print();
    if (options.get_bool("ERI")) mints.ao_eri()->print();

    // return the wavefunction
    return wfn;
}
