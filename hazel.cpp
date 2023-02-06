#include "include/hartreefock.h"
#include "include/molecule.h"
#include <boost/json.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <filesystem>

namespace js = boost::json;
namespace po = boost::program_options;

int main(int argc, char** argv) {
    // initialize the argument parser and container for the arguments
    po::options_description desc("options");
    po::positional_options_description pos;
    po::variables_map vm;

    // add options to the parser
    desc.add_options()
        ("input", po::value<std::string>(), "input file")
        ("nthreads,n", po::value<int>()->default_value(1), "number of threads to use")
        ("version,v", "print version string")
        ("help,h", "produce help message")
    ;pos.add("input", 1);

    // extract the variables from the command line
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm); po::notify(vm);

    // print help if the help flag was provided
    if (vm.count("help")) {
        std::cout << desc << std::endl; return EXIT_SUCCESS;
    }

    // open the provided JSON input
    std::ifstream file(vm["input"].as<std::string>());
    std::stringstream buffer; buffer << file.rdbuf();
    js::value input = js::parse(buffer.str());

    #if defined(_OPENMP)
    libint2::nthreads = vm["nthreads"].as<int>(),
    omp_set_num_threads(libint2::nthreads);
    #endif

    // extract strings from the JSON input
    std::string path = std::filesystem::absolute(vm["input"].as<std::string>()).parent_path();
    std::string molfile = js::value_to<std::string>(input.at("molecule"));
    std::string basis = js::value_to<std::string>(input.at("basis"));
    
    // initialize Hartree-Fock options
    HartreeFockOptions opt = {
        js::value_to<double>(input.at("method").at("thresh")),
        js::value_to<int>(input.at("method").at("maxiter")),
        {
            js::value_to<int>(input.at("method").at("diis").at("start")),
            js::value_to<int>(input.at("method").at("diis").at("keep")),
            js::value_to<bool>(input.at("method").at("diis").at("enabled")),
            js::value_to<double>(input.at("method").at("diis").at("damp"))
        }
    };

    // start the timer
    auto start = Timer::now();

    // print the initial info
    std::cout << "HAZEL" << std::endl << std::endl;

    // initialize the molecule and Hartree-Fock method
    Molecule molecule(path + "/" + molfile, basis);
    HartreeFock hfock(opt);

    // perform the SCF cycle
    libint2::initialize();
    auto result = hfock.scf(molecule);
    libint2::finalize();

    // print orbital energies
    std::cout << "ITER OCC        E [Eh]                E [eV]" << std::endl;
    for (int i = 0; i < result.Eo.rows(); i++) {
        std::cout << boost::format("%4i %3.1f %20.14f %22.14f") % i % ((i + 1) <= result.nocc ? 2.0 : 0.0) % result.Eo(i) % (result.Eo(i) * EH2EV) << "\n";
    }
    std::cout << std::endl;

    // print final energy
    std::cout << boost::format("FINAL SINGLE POINT ENERGY: %.14f Eh") % result.E << std::endl << std::endl;

    // print elapsed time
    std::cout << boost::format("DONE IN %s") % Timer::format(Timer::elapsed(start)) << std::endl;
}
