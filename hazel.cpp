#include "include/hartreefock.h"
#include "include/molecule.h"
#include <boost/json.hpp>
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
        ("ncores,n", po::value<int>()->default_value(1), "number of CPU cores to use")
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
    libint2::nthreads = vm["ncores"].as<int>(),
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
    Printer::printTitle();
    
    // initialize the molecule and Hartree-Fock method
    Molecule molecule(path + "/" + molfile, basis);
    HartreeFock hfock(opt);

    // perform the SCF cycle
    libint2::initialize();
    auto result = hfock.scf(molecule);
    libint2::finalize();

    // print the results and elapsed time
    Printer::printResult(result);
    /* Printer::printTimings(result); */
    Printer::printElapsed(Timer::elapsed(start));
}
