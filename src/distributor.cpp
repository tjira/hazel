#include "distributor.h"
#include <cctype>

// #define TOLOWER(V) [](auto v) {for (auto& e : v) {std::transform(e.begin(), e.end(), e.begin(), [](auto c){return std::tolower(c);});} return v;}(V)
#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // add positional arguments to the argument parser
    program.add_argument("system").help("-- Quantum system to use in .xyz format.").default_value("molecule.xyz");

    // add optional arguments
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-d", "--diis").help("-- Start iteration and Fock history length for DIIS.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    program.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(false).implicit_value(true);
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

    // set the path where to find all the basis functions
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        std::string path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
    }
    
    // set cout global flags
    std::cout << std::fixed << std::setprecision(14);
}

Distributor::~Distributor() {
    std::cout << "\n" + std::string(104, '-') << std::endl;
    std::printf("TOTAL EXECUTION TIME: %s\n", Timer::Format(Timer::Elapsed(start)).c_str());
    std::cout << std::string(104, '-') << std::endl;
}

void Distributor::run() {
    // extract the command line options
    std::pair<int, int> diis = {program.get<std::vector<int>>("-d").at(0), program.get<std::vector<int>>("-d").at(1)};
    int maxiter = program.get<int>("-m"); double thresh = program.get<double>("-t");

    // print the title
    std::cout << "QUANTUM HAZEL" << std::endl;

    // load the system file and extract printing options
    std::vector<std::string> print = program.get<std::vector<std::string>>("-p");
    System system(program.get("system"), program.get("-b"), 0, 1);

    // transform print vector to lowercase
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // print the info header
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n";
    std::cout << std::string(104, '-') + "\n\n";

    // print the info
    std::printf("TIMESTAMP: %s\nCOMPILE FLAGS: %s\n", __TIMESTAMP__, CXXFLAGS);

    // print the system header
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + system.basis + ")\n";
    std::cout << std::string(104, '-') + "\n\n";

    // print the system settings
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);

    // print the system coordinates and distance matrix
    std::cout << "\nSYSTEM COORDINATES\n" << system.coords << std::endl; 
    std::cout << "\nDISTANCE MATRIX\n" << system.dists << std::endl; 

    // define integral struct and initialize libint
    Integrals ints; libint2::initialize();

    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(ints.S = Integral::Overlap(system))
    if (CONTAINS(print, "s")) std::cout << "\n" << ints.S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(ints.T = Integral::Kinetic(system))
    if (CONTAINS(print, "t")) std::cout << "\n" << ints.T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(ints.V = Integral::Nuclear(system))
    if (CONTAINS(print, "v")) std::cout << "\n" << ints.V << std::endl;

    // calculate the electron-electron repulsion integral
    std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(ints.J = Integral::Coulomb(system))
    if (CONTAINS(print, "j")) {std::cout << "\n" << ints.J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.get<bool>("-g") || program.get<bool>("-o")) {
        // calculate the overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(ints.dS = Integral::dOverlap(system))
        if (CONTAINS(print, "ds")) std::cout << "\n" << ints.dS << std::endl;

        // calculate the kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(ints.dT = Integral::dKinetic(system))
        if (CONTAINS(print, "dt")) std::cout << "\n" << ints.dT << std::endl;

        // calculate the nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(ints.dV = Integral::dNuclear(system))
        if (CONTAINS(print, "dv")) std::cout << "\n" << ints.dV << std::endl;

        // calculate the electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(ints.dJ = Integral::dCoulomb(system))
        if (CONTAINS(print, "dj")) {std::cout << "\n" << ints.dJ;} std::cout << "\n";
    }

    // create the initial guess for the density matrix and define the gradient
    Matrix D = Matrix::Zero(ints.S.rows(), ints.S.cols()), G = Matrix::Zero(system.atoms.size(), 3);

    // optimize the molecule
    if (program.get<bool>("-o")) {
        // print the analytical RHF gradient method header
        std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION \n";
        std::cout << std::string(104, '-') + "\n";

        // calculate and print
        std::tie(system, ints, D, G) = Roothaan(system, maxiter, thresh, diis).optimize(ints);

        // print the new system coordinates and distance matrix
        std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
        std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl; 
    }

    // finalize libint
    libint2::finalize();

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n";
    std::cout << std::string(104, '-') + "\n\n";

    // print the RHF options
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", maxiter, thresh, diis.first, diis.second);

    // perform the Hartree-Fock method and extract the results
    auto[C, eps, E] = Roothaan(system, maxiter, thresh, diis).scf(ints, D);

    // print the results
    if (CONTAINS(print, "c")) std::cout << "\nCOEFFICIENT MATRIX\n" << C << std::endl;
    if (CONTAINS(print, "d")) std::cout << "\nDENSITY MATRIX\n" << D << std::endl;
    if (CONTAINS(print, "eps")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(eps) << std::endl;

    // print the energy
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL SINGLE POINT ENERGY: " << E << std::endl;

    // calculate the nuclear gradient
    if (program.get<bool>("-g")) {
        // print the analytical RHF gradient method header
        std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK \n";
        std::cout << std::string(104, '-') + "\n\n";

        // calculate the gradient
        if (!program.get<bool>("-o")) {
            G = Roothaan(system, maxiter, thresh, diis).gradient(ints, C, eps);
        }

        // print the gradient and norm
        std::cout << Matrix(G) << "\n\nGRADIENT NORM: ";
        std::printf("%.2e\n", G.norm());
    }
}
