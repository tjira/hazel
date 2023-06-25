#include "distributor.h"
#include "integral.h"

#define LOCALTIME [](){auto t = std::time(nullptr); auto tm = *std::localtime(&t); std::stringstream ss; ss << std::put_time(&tm, "%a %b %e %T %Y"); return ss.str();}()
#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // initialize subcommands
    hf = argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none); mp2 = argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none);
    ci = argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none);

    // add positional arguments to the main argument parser
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-f", "--file").help("-- Quantum system to use in .xyz file format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads to use.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("-e", "--export").help("-- Export matrices and tensor to files.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("--no-coulomb").help("-- Disable calculation of the coulomb tensor.").default_value(false).implicit_value(true);

    // add positional arguments to the HF argument parser
    hf.add_argument("-d", "--diis").help("-- Start iteration and Fock history length for DIIS.").default_value(std::vector<int>{3, 5}).nargs(1, 2).scan<'i', int>();
    hf.add_argument("-e", "--export").help("-- Export matrices and tensor to files.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
    hf.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    hf.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    hf.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(false).implicit_value(true);
    hf.add_argument("-t", "--thresh").help("-- Threshold for conververgence of iterative calculations.").default_value(1e-8).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    mp2.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    mp2.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    mp2.add_argument("-e", "--export").help("-- Export matrices and tensor to files.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the ci argument parser
    ci.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    ci.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    ci.add_argument("-e", "--export").help("-- Export matrices and tensor to files.").default_value<std::vector<std::string>>({}).append();

    // add the parsers
    program.add_subparser(hf); hf.add_subparser(mp2), hf.add_subparser(ci);

    // parse the arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help if requested
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used(hf) && hf.get<bool>("-h")) {
        std::cout << hf.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used(hf) && hf.is_subcommand_used(mp2) && mp2.get<bool>("-h")) {
        std::cout << mp2.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used(hf) && hf.is_subcommand_used(ci) && ci.get<bool>("-h")) {
        std::cout << ci.help().str(); exit(EXIT_SUCCESS);
    }

    // set the path where to find all the basis functions
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        std::string path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
    }
    
    // set cout global flags and number of threads
    std::cout << std::fixed << std::setprecision(14);
    nthread = program.get<int>("--nthread");
}

Distributor::~Distributor() {
    std::cout << "\n" + std::string(104, '-') << std::endl;
    std::printf("TOTAL EXECUTION TIME: %s\n", Timer::Format(Timer::Elapsed(start)).c_str());
    std::cout << std::string(104, '-') << std::endl;
}


void Distributor::run() {
    // initialize the system
    System system(program.get("-f"), program.get("-b"), 0, 1);

    // extract printing and saving options
    print = program.get<std::vector<std::string>>("-p"), save = program.get<std::vector<std::string>>("-e");

    // transform the print and save vectors to lowercase
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : save) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // print the title
    std::cout << "QUANTUM HAZEL" << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, LOCALTIME.c_str());
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
    std::printf("\nCOMPILER VERSION: GCC %d.%d.%d\nCOMPILER FLAGS: %s\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nTHREADS: %d\n", nthread);

    // print the system block
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + system.basis + ")\n" << std::string(104, '-') + "\n\n";
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);
    std::cout << "\nSYSTEM COORDINATES\n" << system.coords << std::endl; 

    // print the distances if requested
    if (CONTAINS(print, "dist")) std::cout << "\nDISTANCE MATRIX\n" << system.dists << std::endl; 

    // create the initial guess for the density matrix, define the gradient and calculate integrals
    Matrix D = Matrix::Zero(system.shells.nbf(), system.shells.nbf());
    Matrix G = Matrix::Zero(system.atoms.size(), 3), C;
    Tensor<4> JMO; Vector eps; double E, Ecorr;

    // perform the HF calculation
    if (program.is_subcommand_used("hf")) {
        // calculate the molecular integrals
        Integrals ints = integrals(system);

        // extract HF options
        std::pair<int, int> diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
        auto hfprint = hf.get<std::vector<std::string>>("-p"), hfsave = hf.get<std::vector<std::string>>("-e");
        int maxiter = hf.get<int>("-m"); double thresh = hf.get<double>("-t");

        // transform the print and save vectors to lowercase
        for (auto& element : hfprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
        for (auto& element : hfsave) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

        // optimize the molecule with HF method
        if (hf.is_used("-o")) {
            // print the analytical RHF optimization method header and optimize the molecule
            std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION \n" << std::string(104, '-') + "\n";
            std::tie(system, ints, D, G) = Roothaan(system, maxiter, thresh, diis).optimize(ints);

            // print the new system coordinates
            std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 

            // print the new distance matrix
            if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl; 
        }

        // print the RHF method header
        std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";

        // print the RHF options
        std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", maxiter, thresh, diis.first, diis.second);

        // perform the Hartree-Fock method and extract the results
        std::tie(C, eps, E) = Roothaan(system, maxiter, thresh, diis).scf(ints, D);

        // print the results
        if (CONTAINS(hfprint, "eps") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(eps) << std::endl;
        if (CONTAINS(hfprint, "c") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX\n" << C << std::endl;
        if (CONTAINS(hfprint, "d") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX\n" << D << std::endl;

        // save the results
        if (CONTAINS(hfsave, "eps") || CONTAINS(save, "all")) Eigen::Write("EPS.mat", eps);
        if (CONTAINS(hfsave, "c") || CONTAINS(save, "all")) Eigen::Write("C.mat", C);
        if (CONTAINS(hfsave, "d") || CONTAINS(save, "all")) Eigen::Write("D.mat", D);

        // print the energy
        std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
        std::cout << "FINAL HARTREE-FOCK ENERGY: " << E << std::endl;

        // calculate the nuclear gradient
        if (hf.is_used("-g")) {
            // print the analytical RHF gradient method header
            std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK \n";
            std::cout << std::string(104, '-') + "\n\n";

            // calculate the gradient
            if (!hf.is_used("-o")) {
                G = Roothaan(system, maxiter, thresh, diis).gradient(ints, C, eps);
            }

            // print the gradient and norm
            std::cout << Matrix(G) << "\n\nGRADIENT NORM: ";
            std::printf("%.2e\n", G.norm());
        }

        // calculate the MP2 correlation
        if (hf.is_subcommand_used("mp2")) {
            // throw an error if no coulomb
            if (program.get<bool>("--no-coulomb")) throw std::runtime_error("I'm sorry, you need the coulomb tensor for the MP2 method.");

            // get the printing options
            auto mp2print = mp2.get<std::vector<std::string>>("-p"), mp2save = mp2.get<std::vector<std::string>>("-e");

            // transform the print and save vectors to lowercase
            for (auto& element : mp2print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
            for (auto& element : mp2save) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

            // print the analytical MP2 correlation method header
            std::cout << "\n" + std::string(104, '-') + "\nMP2 CORRELATION ENERGY\n";
            std::cout << std::string(104, '-') + "\n";

            // transform the coulomb tensor
            std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(JMO = Transform::Coulomb(ints.J, C))
            if (CONTAINS(mp2print, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << JMO;} std::cout << "\n";
            if (CONTAINS(mp2save, "jmo") || CONTAINS(save, "all")) Eigen::Write("JMO.mat", JMO);

            // calculate the correlation energy
            Ecorr = MP(system).mp2(JMO, eps);

            // print the gradient and norm
            std::cout << "\nMP2 CORRELATION ENERGY: " << Ecorr << std::endl;
            std::cout << "FINAL MP2 ENERGY: " << E + Ecorr << std::endl;

        // calculate ci correlation
        } else if (hf.is_subcommand_used("ci")) {
            // throw an error if no coulomb
            if (program.get<bool>("--no-coulomb")) throw std::runtime_error("I'm sorry, you need the coulomb tensor for the ccsd method.");

            // get the printing options
            auto ciprint = ci.get<std::vector<std::string>>("-p"), cisave = ci.get<std::vector<std::string>>("-e");

            // transform the print and save vectors to lowercase
            for (auto& element : ciprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
            for (auto& element : cisave) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

            // print the CI correlation method header
            std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n";
            std::cout << std::string(104, '-') + "\n";

            // transform the coulomb tensor
            std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(JMO = Transform::Coulomb(ints.J, C))
            if (CONTAINS(ciprint, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << JMO;} std::cout << "\n";
            if (CONTAINS(cisave, "jmo") || CONTAINS(save, "all")) Eigen::Write("JMO.mat", JMO);

            // calculate the correlation energy retrieved
            auto[CIH, CIC, CIE] = CI(system, E - Integral::Repulsion(system)).cid(ints, JMO, C);
            double Ecorr = CIE(0) - E + Integral::Repulsion(system);

            // print the result matrices
            if (CONTAINS(ciprint, "cih") || CONTAINS(print, "all")) {std::cout << "\n" << CIH << "\n";}
            if (CONTAINS(ciprint, "cic") || CONTAINS(print, "all")) {std::cout << "\n" << CIC << "\n";}
            if (CONTAINS(ciprint, "cie") || CONTAINS(print, "all")) {std::cout << "\n" << CIE << "\n";}

            // save the result matrices
            if (CONTAINS(cisave, "cih") || CONTAINS(save, "all")) Eigen::Write("CIH.mat", CIH);
            if (CONTAINS(cisave, "cic") || CONTAINS(save, "all")) Eigen::Write("CIC.mat", CIC);
            if (CONTAINS(cisave, "cie") || CONTAINS(save, "all")) Eigen::Write("CIE.mat", CIE);

            // print the gradient and norm
            std::cout << "\nCI CORRELATION ENERGY: " << Ecorr << std::endl;
            std::cout << "FINAL CI ENERGY: " << E + Ecorr << std::endl;
        }
    }
}

Integrals Distributor::integrals(const System& system) const {
    // define the struct
    Integrals ints;

    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(ints.S = Integral::Overlap(system))
    if (CONTAINS(print, "s") || CONTAINS(print, "all")) std::cout << "\n" << ints.S << std::endl;
    if (CONTAINS(save, "s") || CONTAINS(save, "all")) Eigen::Write("S.mat", ints.S);

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(ints.T = Integral::Kinetic(system))
    if (CONTAINS(print, "t") || CONTAINS(print, "all")) std::cout << "\n" << ints.T << std::endl;
    if (CONTAINS(save, "t") || CONTAINS(save, "all")) Eigen::Write("T.mat", ints.T);

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(ints.V = Integral::Nuclear(system))
    if (CONTAINS(print, "v") || CONTAINS(print, "all")) std::cout << "\n" << ints.V << std::endl;
    if (CONTAINS(save, "v") || CONTAINS(save, "all")) Eigen::Write("V.mat", ints.V);

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(ints.J = Integral::Coulomb(system))}
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(print, "j") || CONTAINS(print, "all"))) {std::cout << "\n" << ints.J;} std::cout << "\n";
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(save, "j") || CONTAINS(save, "all"))) Eigen::Write("J.mat", ints.J);

    // if derivatives of the integrals are needed
    if (hf.is_used("-g") || hf.is_used("-o")) {
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

    // return integrals
    return ints;
}
