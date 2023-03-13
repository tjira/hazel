#include "include/hartreefock.h"
#include "include/defaults.h"
#include <argparse/argparse.hpp>

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options::PRINT, kinetic, oneelec, overlap, density, orben, mos, iter);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options::MDYN, timestep, steps, output);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options::GRAD::PRINT, iter);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options::GRAD, increment, nthread, print);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options::DIIS, start, keep, damp);
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(HartreeFock::Options, thresh, maxiter, diis, engrad, print, mulliken);

std::tuple<json, argparse::ArgumentParser> parse(int argc, char** argv) {
    // initialize the argument parser and container for the arguments
    argparse::ArgumentParser program("Hazel", "1.0", argparse::default_arguments::none);

    // add options to the parser
    program.add_argument("input").help("Hazel input file.");
    program.add_argument("-h").help("Display this help message and exit.").default_value(false).implicit_value(true);

    // extract the variables from the command line
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl << std::endl << program; exit(EXIT_FAILURE);
    }

    // print help if the help flag was provided
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    }

    // check if the input file exists
    if (!std::filesystem::exists(program.get<std::string>("input"))) {
        throw std::runtime_error("Input file does not exist.");
    }

    // create the input object
    json input = json::parse(std::ifstream(program.get<std::string>("input")));

    // return the input
    return { input, program };
}

json patch(json input) {
    if (input.at("method").at("name") == "HF") Defaults::hfopt.merge_patch(input.at("method")), input.at("method") = Defaults::hfopt;
    if (input.contains("dynamics")) Defaults::mdopt.merge_patch(input.at("dynamics")), input.at("dynamics") = Defaults::mdopt;
    if (input.contains("print")) { Defaults::print.merge_patch(input.at("print")); } input["print"] = Defaults::print;
    return input;
}

int main(int argc, char** argv) {
    // parse the command line arguments and create a logger
    auto [input, program] = parse(argc, argv);

    // print the initial info
    Logger::Log(false, "HAZEL\nCOMPILE FLAGS: %s", STRINGIFY(GPPFLAGS));

    // print the input if requested
    Logger::Log(!patch(input).at("print").at("input"), "\nINPUT\n%s", input.dump(4));

    // initialize the system
    std::string basis = input.at("basis").get<std::string>();
    System system(PABSOLUTE(input.at("system").get<std::string>()), basis);

    // print the system
    Logger::Log(!patch(input).at("print").at("system"), "\nSYSTEM FILE\n%s", FCONTENTS(PABSOLUTE(input.at("system").get<std::string>())));

    // print the system specification
    Logger::Log(!patch(input).at("print").at("molecule"), "\nMOLECULE\nATOMS: %i, ELECTRONS: %i, BASIS: %i", system.getSize(), system.getElectrons(), basis);

    // initialize the libint library nd start the timer
    libint2::initialize(); auto start = Timer::now();

    // initialize method options
    HartreeFock::Options opt = patch(input).at("method").get<HartreeFock::Options>();
    opt.diis.on = input.at("method").contains("diis");

    // insert dynamics struct if requested
    if (input.contains("dynamics")) {
        opt.dyn = patch(input).at("dynamics").get<HartreeFock::Options::MDYN>();
        opt.dyn.output = PABSOLUTE(opt.dyn.output);
    }

    // initialize HF and perform the SCF cycle
    HartreeFock hfock(opt); auto result = hfock.scf(system, false);

    // calculate and print the energy gradient if requested
    if (input.at("method").contains("engrad") || input.contains("dynamics")) HartreeFock(opt).gradient(system, false);

    // perform the mulliken analysis if requested
    if (patch(input).at("method").at("mulliken")) {
        auto mulliken = system.mulliken(result.D);
        Logger::Log(false, "\nMULLIKEN CHARGES\n%=2s %=9s", "SM", "Q");
        for (size_t i = 0; i < system.getAtoms().size(); i++) {
            Logger::Log(false, "%2s %9.6f", an2sm.at(system.getAtoms().at(i).atomic_number), mulliken.q(i));
        }
        Logger::Log(false, "SUM OF MULLIKEN CHARGES: %.6f", mulliken.q.sum());
    }

    // perform the molecular dynamics if requested
    if (input.contains("dynamics")) hfock.dynamics(system, false);

    // print elapsed time
    Logger::Log(!patch(input).at("print").at("time"), "\nDONE IN %s", Timer::format(Timer::elapsed(start)));

    // finalize the libint library 
    libint2::finalize();
}
