#include "include/hartreefock.h"
#include "include/moleculardynamics.h"
#include <boost/json/src.hpp>
#include <boost/format.hpp>
#include <boost/json/value_to.hpp>
#include <boost/program_options.hpp>
#include <filesystem>

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X

namespace po = boost::program_options;
namespace js = boost::json;

js::value fill(js::value value) {
    if (!value.at("method").as_object().if_contains("diis")) {
        value.at("method").as_object().insert(std::pair<std::string, js::object>("diis", {}));
    }
    value.at("method").as_object().insert(std::pair<std::string, const char*>("output-trajectory", "trajectory.xyz"));
    value.at("method").as_object().insert(std::pair<std::string, int>("maxiter", 100));
    value.at("method").as_object().insert(std::pair<std::string, double>("thresh", 1e-8));
    value.at("method").at("diis").as_object().insert(std::pair<std::string, bool>("enabled", true));
    value.at("method").at("diis").as_object().insert(std::pair<std::string, double>("damp", 0.0));
    value.at("method").at("diis").as_object().insert(std::pair<std::string, int>("start", 3));
    value.at("method").at("diis").as_object().insert(std::pair<std::string, int>("keep", 5));
    return value;
}

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
    if (!std::filesystem::exists(vm["input"].as<std::string>())) {
        throw std::runtime_error("Input file does not exist.");
    }
    std::ifstream file(vm["input"].as<std::string>());
    std::stringstream buffer; buffer << file.rdbuf();
    js::value input = fill(js::parse(buffer.str()));

    // print the initial info and input
    std::cout << "HAZEL\nCOMPILE FLAGS: " << STRINGIFY(GPPFLAGS) << "\n" << std::endl;
    std::cout << "INPUT\n" << buffer.str() << std::endl;

    // set number of threads
    libint2::nthreads = vm["nthreads"].as<int>();
    #if defined(_OPENMP)
    omp_set_num_threads(libint2::nthreads);
    #endif

    // extract strings from the JSON input
    std::string path = std::filesystem::absolute(vm["input"].as<std::string>()).parent_path();
    std::string sysfile = js::value_to<std::string>(input.at("system"));

    // if Hartree-Fock
    if (js::value_to<std::string>(input.at("method").at("name")) == "HF") {
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

        // initialize the libint library 
        libint2::initialize();

        // initialize the system
        if (!std::filesystem::exists(path + "/" + sysfile)) {
            throw std::runtime_error("System file does not exist.");
        }
        System system(path + "/" + sysfile, basis);

        // print the system specification
        std::cout << boost::format("MOLECULE\nATOMS: %i, ELECTRONS: %i, BASIS: %i\n") % system.getAtomCount() % system.getElectronCount() % basis << std::endl;

        // initialize HF
        HartreeFock hfock(opt);

        // perform the SCF cycle
        auto result = hfock.scf(system);

        // perform the mulliken analysis
        auto mulliken = system.mulliken(result.D);

        // print orbital energies
        std::cout << "ORBITAL ENERGIES AND OCCUPATION\nITER OCC        E [Eh]                E [eV]\n";
        for (int i = 0; i < result.eps.rows(); i++) {
            std::cout << boost::format("%4i %3.1f %20.14f %22.14f\n") % i % ((i + 1) <= result.nocc ? 2.0 : 0.0) % result.eps(i) % (result.eps(i) * EH2EV);
        }
        std::cout << std::endl;

        // print the mulliken charges
        std::cout << "MULLIKEN CHARGES\nAN     Q\n";
        for (int i = 0; i < system.getAtomCount(); i++) {
            std::cout << boost::format("%2i %9.6f\n") % system.getAtom(i).atomic_number % mulliken.q(i);
        }
        std::cout << "\nSUM OF MULLIKEN CHARGES: " << boost::format("%.6f\n") % mulliken.q.sum() << std::endl;

        // print final energy
        std::cout << boost::format("FINAL SINGLE POINT ENERGY: %.14f Eh\n") % result.E << std::endl;

        // print elapsed time
        std::cout << boost::format("DONE IN %s") % Timer::format(Timer::elapsed(start)) << std::endl;

        // finalize the libint library 
        libint2::finalize();
    } else if (js::value_to<std::string>(input.at("method").at("name")) == "MD") {
        MolecularDynamicsOptions opt {
            js::value_to<double>(input.at("method").at("timestep")),
            js::value_to<int>(input.at("method").at("steps"))
        };

        // start the timer
        auto start = Timer::now();

        // initialize the system
        if (!std::filesystem::exists(path + "/" + sysfile)) {
            throw std::runtime_error("System file does not exist.");
        }
        System system(path + "/" + sysfile);

        MolecularDynamics mdyn(opt);

        MolecularDynamicsResult result = mdyn.run(system.getParticles(), path + "/" + js::value_to<std::string>(input.at("method").at("output-trajectory")));

        // print elapsed time
        std::cout << boost::format("DONE IN %s") % Timer::format(Timer::elapsed(start)) << std::endl;
    }
}
