#include "distributor.h"

#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // add positional arguments to the argument parser
    program.add_argument("system").help("-- Quantum system to use in .xyz format.").default_value("molecule.xyz");

    // add optional arguments
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-d", "--diis").help("-- Start iteration and Fock history length for DIIS.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    program.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
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
}

Distributor::~Distributor() {
    std::cout << "\n" + std::string(104, '-') << std::endl;
    std::cout << std::format("TOTAL EXECUTION TIME: {:s}", Timer::Format(Timer::Elapsed(start))) << std::endl;
    std::cout << std::string(104, '-') << std::endl;
}

void Distributor::run() {
    // extract the command line options
    std::pair<int, int> diis = {program.get<std::vector<int>>("-d").at(0), program.get<std::vector<int>>("-d").at(1)};
    int maxiter = program.get<int>("-m"); double thresh = program.get<double>("-t");

    // print the title
    std::cout << "HAZEL" << std::endl;

    // load the system file and extract printing options
    std::vector<std::string> print = program.get<std::vector<std::string>>("-p");
    System system(program.get("system"), program.get("-b"), 0, 1);

    // define the integral matrices and tensors
    Tensor<3> dS, dT, dV; Tensor<5> dJ;
    Matrix S, T, V; Tensor<4> J;

    // initialize libint
    libint2::initialize();

    // calculate the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(S = Integral::Overlap(system))
    if (CONTAINS(print, "S")) std::cout << "\n" << S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(T = Integral::Kinetic(system))
    if (CONTAINS(print, "T")) std::cout << "\n" << T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(V = Integral::Nuclear(system))
    if (CONTAINS(print, "V")) std::cout << "\n" << V << std::endl;

    // calculate the electron-electron repulsion integral
    std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(J = Integral::Coulomb(system))
    if (CONTAINS(print, "J")) {std::cout << "\n" << J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.get<bool>("-g")) {
        // calculate the overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(dS = Integral::dOverlap(system))
        if (CONTAINS(print, "dS")) std::cout << "\n" << dS << std::endl;

        // calculate the kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(dT = Integral::dKinetic(system))
        if (CONTAINS(print, "dT")) std::cout << "\n" << dT << std::endl;

        // calculate the nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(dV = Integral::dNuclear(system))
        if (CONTAINS(print, "dV")) std::cout << "\n" << dV << std::endl;

        // calculate the electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(dJ = Integral::dCoulomb(system))
        if (CONTAINS(print, "dJ")) {std::cout << "\n" << dJ;} std::cout << "\n";
    }

    // finalize libint
    libint2::finalize();

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n";
    std::cout << std::format("MAXITER: {:d}, THRESH: {:.2e}, DIIS: [START: {:d}, KEEP: {:d}]\n", maxiter, thresh, diis.first, diis.second);
    std::cout << std::string(104, '-') + "\n\n";

    // create the initial guess for the density matrix
    Matrix D = Matrix::Zero(S.rows(), S.cols());

    // perform the Hartree-Fock method and extract the results
    auto[C, eps, Eel] = Roothaan(system, maxiter, thresh, diis).scf(T + V, J, S, D);

    // print the results
    if (CONTAINS(print, "C")) std::cout << "\nCOEFFICIENT MATRIX\n" << C << std::endl;
    if (CONTAINS(print, "D")) std::cout << "\nDENSITY MATRIX\n" << D << std::endl;
    if (CONTAINS(print, "EPS")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(eps) << std::endl;

    // print the energy
    std::cout << std::format("\nTOTAL NUCLEAR REPULSION ENERGY: {:.14f}", Integral::Repulsion(system)) << std::endl;
    std::cout << std::format("FINAL SINGLE POINT ENERGY: {:.14f}", Eel + Integral::Repulsion(system)) << std::endl;

    // calculate the nuclear gradient
    if (program.get<bool>("-g")) {
        // print the analytical RHF gradient method header
        std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK \n";
        std::cout << std::string(104, '-') + "\n\n";

        // calculate and print
        Matrix G = Roothaan(system, maxiter, thresh, diis).gradient(dT, dV, dJ, dS, C, eps);
        std::cout << Matrix(G + Integral::dRepulsion(system)) << std::endl;
    }

}
