#include "include/hartreefock.h"
#include "include/moleculardynamics.h"
#include <argparse/argparse.hpp>
#include <filesystem>

#define STRINGIFY(X) STRING(X)
#define STRING(X) #X

json patch(json input) {
    std::string method = input.at("method").at("name");
    if (method == "HF") Defaults::hfopt.merge_patch(input.at("method")), input.at("method") = Defaults::hfopt;
    if (method == "MD") Defaults::mdopt.merge_patch(input.at("method")), input.at("method") = Defaults::mdopt;
    return input;
}

int main(int argc, char** argv) {
    // initialize the argument parser and container for the arguments
    argparse::ArgumentParser program("Hazel", "1.0", argparse::default_arguments::none);

    // add options to the parser
    program.add_argument("input").help("Hazel input file.");
    program.add_argument("-h").help("Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-n").help("Number of threads to use.").default_value(1).scan<'i', int>();

    // extract the variables from the command line
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl << std::endl << program; return EXIT_FAILURE;
    }

    // print help if the help flag was provided
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); return EXIT_SUCCESS;
    }

    // open the provided JSON input
    if (!std::filesystem::exists(program.get<std::string>("input"))) {
        throw std::runtime_error("Input file does not exist.");
    }
    json input = patch(json::parse(std::ifstream(program.get<std::string>("input"))));

    // print the initial info and input
    std::cout << "HAZEL\nCOMPILE FLAGS: " << STRINGIFY(GPPFLAGS) << "\n" << std::endl;
    std::cout << "INPUT\n" << input.dump(4) << "\n" << std::endl;

    // set number of threads
    libint2::nthreads = program.get<int>("n");
    #if defined(_OPENMP)
    omp_set_num_threads(libint2::nthreads);
    #endif

    // extract strings from the JSON input
    std::string path = std::filesystem::absolute(program.get<std::string>("input")).parent_path();
    std::string sysfile = input.at("system").get<std::string>();

    // if Hartree-Fock
    if (input.at("method").at("name").get<std::string>() == "HF") {
        // initialize basis set string
        std::string basis = input.at("basis").get<std::string>();

        // initialize Hartree-Fock options
        HartreeFockOptions opt = input.at("method").get<HartreeFockOptions>();
        opt.diis = input.at("method").at("diis").get<HartreeFockOptions::DIIS>();

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

    } else if (input.at("method").at("name").get<std::string>() == "MD") {
        // initialize MD options
        MolecularDynamicsOptions opt = input.at("method").get<MolecularDynamicsOptions>();

        // start the timer
        auto start = Timer::now();

        // initialize the system
        if (!std::filesystem::exists(path + "/" + sysfile)) {
            throw std::runtime_error("System file does not exist.");
        }
        System system(path + "/" + sysfile);
        
        // initialize the pair potential
        auto pairOpt = input.at("method").at("potential").at("pair").get<PotentialOptions>();
        auto bondOpt = input.at("method").at("potential").at("bond").get<PotentialOptions>();

        // initialize the molecular dynamics object
        MolecularDynamics mdyn(ForceField(pairOpt, bondOpt, system), opt);

        // perform the molecular dynamics
        MolecularDynamicsResult result = mdyn.run(system.getParticles(), path + "/" + input.at("method").at("output-trajectory").get<std::string>());

        // print elapsed time
        std::cout << boost::format("DONE IN %s") % Timer::format(Timer::elapsed(start)) << std::endl;
    }
}
