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

    // transform print vector to lowercase and create auxiliary system matrices (coordinate and distance matrix)
    for (auto& el : print) std::transform(el.begin(), el.end(), el.begin(), [](auto c){return std::tolower(c);});
    Matrix coords(system.atoms.size(), 3), dists(system.atoms.size(), system.atoms.size());

    // fill the coordinate and distance matrices
    for (int i = 0; i < coords.rows(); i++) coords.row(i) = BOHR2A * Eigen::Vector<double, 3>(system.atoms.at(i).x, system.atoms.at(i).y, system.atoms.at(i).z);
    for (int i = 0; i < coords.rows(); i++) for (int j = 0; j < coords.rows(); j++) dists(i, j) = (coords.row(i) - coords.row(j)).norm();

    // print the info header
    std::cout << "\n" + std::string(104, '-') + "\n";
    std::printf("TIMESTAMP: %s\nCOMPILE FLAGS: %s\n", __TIMESTAMP__, CXXFLAGS);
    std::cout << std::string(104, '-') + "\n";

    // print the system header
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION\n";
    std::printf("ATOMS: %d, CHARGE: %d, MULTIPLICITY: %d, ELECTRONS: %d, BASIS: %s, NBF: %d\n", (int)system.atoms.size(), system.charge, system.multi, system.electrons, program.get("-b").c_str(), (int)system.shells.nbf());
    std::cout << std::string(104, '-') + "\n";

    // print the system coordinates and distance matrix
    std::cout << "\nSYSTEM COORDINATES\n" << coords << std::endl; 
    std::cout << "\nDISTANCE MATRIX\n" << dists << std::endl; 

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
    if (CONTAINS(print, "s")) std::cout << "\n" << S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(T = Integral::Kinetic(system))
    if (CONTAINS(print, "t")) std::cout << "\n" << T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(V = Integral::Nuclear(system))
    if (CONTAINS(print, "v")) std::cout << "\n" << V << std::endl;

    // calculate the electron-electron repulsion integral
    std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(J = Integral::Coulomb(system))
    if (CONTAINS(print, "j")) {std::cout << "\n" << J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.get<bool>("-g")) {
        // calculate the overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(dS = Integral::dOverlap(system))
        if (CONTAINS(print, "ds")) std::cout << "\n" << dS << std::endl;

        // calculate the kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(dT = Integral::dKinetic(system))
        if (CONTAINS(print, "dt")) std::cout << "\n" << dT << std::endl;

        // calculate the nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(dV = Integral::dNuclear(system))
        if (CONTAINS(print, "dv")) std::cout << "\n" << dV << std::endl;

        // calculate the electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(dJ = Integral::dCoulomb(system))
        if (CONTAINS(print, "dj")) {std::cout << "\n" << dJ;} std::cout << "\n";
    }

    // finalize libint
    libint2::finalize();

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n";
    std::printf("MAXITER: %d, THRESH: %.2e, DIIS: [START: %d, KEEP: %d]\n", maxiter, thresh, diis.first, diis.second);
    std::cout << std::string(104, '-') + "\n\n";

    // create the initial guess for the density matrix
    Matrix D = Matrix::Zero(S.rows(), S.cols());

    // perform the Hartree-Fock method and extract the results
    auto[C, eps, Eel] = Roothaan(system, maxiter, thresh, diis).scf(T + V, J, S, D);

    // print the results
    if (CONTAINS(print, "c")) std::cout << "\nCOEFFICIENT MATRIX\n" << C << std::endl;
    if (CONTAINS(print, "d")) std::cout << "\nDENSITY MATRIX\n" << D << std::endl;
    if (CONTAINS(print, "eps")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(eps) << std::endl;

    // print the energy
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL SINGLE POINT ENERGY: " << Eel + Integral::Repulsion(system) << std::endl;

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
