#include "integral.h"
#include "roothaan.h"
#include "argparse.hpp"
#include "timer.h"

#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))

argparse::ArgumentParser parse(int argc, char** argv) {
    // create the argument parser
    argparse::ArgumentParser program("hazel", "0.1", argparse::default_arguments::none);

    // add positional arguments
    program.add_argument("system").help("-- Quantum system to use in .xyz format.").default_value("molecule.xyz");

    // arr optional arguments
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("-t", "--thresh").help("-- Threshold for conververgence of iterative calculations.").default_value(1e-12).scan<'g', double>();

    // parse the arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help if requested
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    }

    //return the program
    return program;
}

int main(int argc, char** argv) {
    // parse the arguments
    argparse::ArgumentParser program = parse(argc, argv);

    // print the title
    std::cout << "HAZEL" << std::endl;

    // load the system file and extract printing options
    std::vector<std::string> print = program.get<std::vector<std::string>>("-p");
    System system(program.get("system"), program.get("-b"), 0, 1);

    // initialize libint and start the timer
    libint2::initialize(); Timer::Timepoint start = Timer::Now();

    // calculate the integral calculation header
    if (CONTAINS(print, "S") || CONTAINS(print, "T") || CONTAINS(print, "V")) std::cout << "\n" + std::string(WIDTH, '-') + "\nINTEGRAL CALCULATION\n";
    if (CONTAINS(print, "S") || CONTAINS(print, "T") || CONTAINS(print, "V")) std::cout << std::string(WIDTH, '-') << std::endl;

    // calculate the overlap integral
    if (CONTAINS(print, "S")) std::cout << "\nOVERLAP INTEGRAL: " << std::flush;
    Matrix S = Integral::Overlap(system);
    if (CONTAINS(print, "S")) std::cout << "\n" << S << std::endl;

    // calculate the kinetic integral
    if (CONTAINS(print, "T")) std::cout << "\nKINETIC INTEGRAL: " << std::flush;
    Matrix T = Integral::Kinetic(system);
    if (CONTAINS(print, "T")) std::cout << "\n" << T << std::endl;

    // calculate the nuclear-electron attraction integral
    if (CONTAINS(print, "V")) std::cout << "\nNUCLEAR INTEGRAL: " << std::flush;
    Matrix V = Integral::Nuclear(system);
    if (CONTAINS(print, "V")) std::cout << "\n" << V << std::endl;

    // calculate the electron-electron repulsion integral
    if (CONTAINS(print, "J")) std::cout << "\nCOULOMB INTEGRAL: " << std::flush;
    Tensor<4> J = Integral::Coulomb(system); 
    if (CONTAINS(print, "J")) std::cout << "\n" << J << std::endl;

    // perform the Hartree-Fock method and extract the results
    auto[C, eps, Eel] = Roothaan(program.get<int>("-m"), program.get<double>("-t")).scf(T + V, J, S, Matrix(S.rows()), system.electrons / 2);

    // print the final energy and execution time
    std::cout << "\n" + std::string(WIDTH, '-') << std::endl;
    std::cout << std::format("TOTAL NUCLEAR REPULSION ENERGY: {:.14f}", Integral::Repulsion(system)) << std::endl;
    std::cout << std::format("FINAL SINGLE POINT ENERGY: {:.14f}", Eel + Integral::Repulsion(system)) << std::endl;
    std::cout << std::format("TOTAL EXECUTION TIME: {:s}", Timer::Format(Timer::Elapsed(start))) << std::endl;
    std::cout << std::string(WIDTH, '-') << std::endl;

    // finalize libint
    libint2::finalize();
}
