#include "include/hartreefock.h"
#include "include/molecule.h"

#include "lib/argparse/argparse.hpp"
#include "lib/json/json.hpp"

int main(int argc, char** argv) {
    argparse::ArgumentParser program("Hazel", "0.1", argparse::default_arguments::none);
    program.add_argument("input").help("Input file in JSON format.");
    Printer::printTitle();

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl << std::endl << program; return EXIT_FAILURE;
    }

    nlohmann::json input = nlohmann::json::parse(std::ifstream(program.get<std::string>("input")));

    HartreeFockOptions opt = { .damp = input.at("damp"), .thresh = input.at("thresh"), .maxiter = input.at("maxiter") };
    HartreeFock hfock(opt);
    Molecule molecule(input.at("molecule"), input.at("basis"));

    libint2::initialize();
    auto result = hfock.scf(molecule);
    libint2::finalize();
}
