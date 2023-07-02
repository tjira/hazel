#include "distributor.h"

#define LOCALTIME [](){auto t = std::time(nullptr); auto tm = *std::localtime(&t); std::stringstream ss; ss << std::put_time(&tm, "%a %b %e %T %Y"); return ss.str();}()
#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // initialize subcommands
    hf = argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none); mp2 = argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none);
    ci = argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none);

    // add positional arguments to the main argument parser
    program.add_argument("-c", "--center").help("-- Center the molecule before doing any calculation.").default_value(false).implicit_value(true);
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-f", "--file").help("-- Quantum system to use in .xyz file format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads to use.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("-e", "--export").help("-- Export matrices and tensors to files.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("--no-coulomb").help("-- Disable calculation of the coulomb tensor.").default_value(false).implicit_value(true);

    // add positional arguments to the HF argument parser
    hf.add_argument("-d", "--diis").help("-- Start iteration and Fock history length for DIIS.").default_value(std::vector<int>{3, 5}).nargs(1, 2).scan<'i', int>();
    hf.add_argument("-e", "--export").help("-- Export matrices and tensors to files.").default_value<std::vector<std::string>>({}).append();
    // hf.add_argument("-f", "--frequency").help("-- Enable frequency calculation.").default_value(false).implicit_value(true);
    hf.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
    hf.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    hf.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    hf.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    hf.add_argument("-t", "--thresh").help("-- Threshold for conververgence.").default_value(1e-12).scan<'g', double>();
    hf.add_argument("--numgrad").help("-- Perform the gradient calculation numerically.").default_value(1e-5).scan<'g', double>();
    // hf.add_argument("--numhess").help("-- Perform the frequency calculation numerically.").default_value(1e-5).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    mp2.add_argument("-e", "--export").help("-- Export matrices and tensors to files.").default_value<std::vector<std::string>>({}).append();
    // mp2.add_argument("-f", "--frequency").help("-- Enable frequency calculation.").default_value(false).implicit_value(true);
    mp2.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
    mp2.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    mp2.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    mp2.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    mp2.add_argument("--numgrad").help("-- Perform the gradient calculation numerically.").default_value(1e-5).scan<'g', double>();
    // mp2.add_argument("--numhess").help("-- Perform the frequency calculation numerically.").default_value(1e-5).scan<'g', double>();

    // add positional arguments to the CI argument parser
    ci.add_argument("-e", "--export").help("-- Export matrices and tensors to files.").default_value<std::vector<std::string>>({}).append();
    // ci.add_argument("-f", "--frequency").help("-- Enable frequency calculation.").default_value(false).implicit_value(true);
    // ci.add_argument("-g", "--gradient").help("-- Enable gradient calculation.").default_value(false).implicit_value(true);
    ci.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    ci.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    // ci.add_argument("--numgrad").help("-- Perform the gradient calculation numerically.").default_value(1e-5).scan<'g', double>();
    // ci.add_argument("--numhess").help("-- Perform the frequency calculation numerically.").default_value(1e-5).scan<'g', double>();

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
    Data data; data.system = System(program.get("-f"), program.get("-b"), 0, 1);

    // extract printing and saving options
    print = program.get<std::vector<std::string>>("-p"), save = program.get<std::vector<std::string>>("-e");
    std::vector<std::string> ciprint, hfprint, mp2print, cisave, mp2save, hfsave;

    // transform the print and save vectors to lowercase
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : save) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // print the title
    std::cout << "QUANTUM HAZEL (" << nthread << " THREAD" << (nthread > 1 ? "S)" : ")") << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, LOCALTIME.c_str());
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
    std::printf("\nCOMPILER VERSION: GCC %d.%d.%d\nCOMPILER FLAGS: %s\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);

    // print the system block
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + data.system.basis + ")\n" << std::string(104, '-') + "\n\n";
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)data.system.atoms.size(), data.system.electrons, (int)data.system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", data.system.charge, data.system.multi);
    std::cout << "\nSYSTEM COORDINATES\n" << data.system.coords << std::endl; 

    // center the molecule if requested
    if (program.get<bool>("--center")) {
        Matrix dir(data.system.atoms.size(), 3); dir.rowwise() -= data.system.coords.colwise().sum() / data.system.atoms.size();
        data.system.move(dir * A2BOHR); std::cout << "\nCENTERED SYSTEM COORDINATES\n" << data.system.coords << std::endl; 
    }

    // print the distances if requested
    if (CONTAINS(print, "dist")) std::cout << "\nDISTANCE MATRIX\n" << data.system.dists << std::endl; 

    // extract all the options
    if (program.is_subcommand_used("hf")) {
        hfprint = hf.get<std::vector<std::string>>("-p"), hfsave = hf.get<std::vector<std::string>>("-e");
        data.roothaan.diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
        data.roothaan.maxiter = hf.get<int>("-m"), data.roothaan.grad.step = hf.get<double>("--numgrad");
        data.roothaan.thresh = hf.get<double>("-t"), data.roothaan.opt.thresh = hf.get<double>("-o");
        data.roothaan.grad.numerical = hf.is_used("--numgrad");
        if (hf.is_subcommand_used("mp2")) {
            mp2print = mp2.get<std::vector<std::string>>("-p"), mp2save = mp2.get<std::vector<std::string>>("-e");
            data.mp.grad.numerical = mp2.is_used("--numgrad"), data.mp.grad.step = mp2.get<double>("--numgrad");
            data.mp.opt.thresh = mp2.get<double>("-o");
        } else if (hf.is_subcommand_used("ci")) {
            ciprint = ci.get<std::vector<std::string>>("-p"), cisave = ci.get<std::vector<std::string>>("-e");
        }
    }

    // transform the print and save vectors to lowercase
    for (auto& element : mp2print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : hfprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : mp2save) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : ciprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : hfsave) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : cisave) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // perform the HF calculation
    if (program.is_subcommand_used("hf")) {
        // calculate the molecular integrals and create the guess density
        data = integrals(data), data.roothaan.D = Matrix::Zero(data.system.shells.nbf(), data.system.shells.nbf());

        // optimize the molecule with HF method
        if (hf.is_used("-o")) {
            // print the analytical RHF optimization method, optimize the molecule and print the new coordinates
            std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n";
            
            // perform the optimization
            data = Roothaan(data).optimize();

            // print the optimized coordinates
            std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << data.system.coords << std::endl; 

            // print the new distance matrix if requested
            if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << data.system.dists << std::endl; 
        }

        // print the RHF method header
        std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";

        // print the RHF options
        std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", data.roothaan.maxiter, data.roothaan.thresh, data.roothaan.diis.start, data.roothaan.diis.keep);

        // perform the Hartree-Fock calculation
        data = Roothaan(data).scf();

        // print and save the resulting matrices
        if (CONTAINS(hfprint, "eps") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(data.roothaan.eps) << std::endl;
        if (CONTAINS(hfprint, "c") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX\n" << data.roothaan.C << std::endl;
        if (CONTAINS(hfprint, "d") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX\n" << data.roothaan.D << std::endl;
        if (CONTAINS(hfsave, "eps") || CONTAINS(save, "all")) Eigen::Write("EPS.mat", data.roothaan.eps);
        if (CONTAINS(hfsave, "c") || CONTAINS(save, "all")) Eigen::Write("C.mat", data.roothaan.C);
        if (CONTAINS(hfsave, "d") || CONTAINS(save, "all")) Eigen::Write("D.mat", data.roothaan.D);

        // print the energy
        std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(data.system) << std::endl;
        std::cout << "FINAL HARTREE-FOCK ENERGY: " << data.roothaan.E << std::endl;

        // calculate the nuclear gradient
        if (hf.is_used("-g")) {
            // print the analytical RHF gradient method header
            if (hf.is_used("--numgrad")) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
            else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
            std::cout << std::string(104, '-') + "\n\n";

            // perform the gradient calculation
            if (!hf.is_used("-o")) data = Roothaan(data).gradient();

            // print the gradient and norm
            std::cout << data.roothaan.grad.G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", data.roothaan.grad.G.norm());
        }

        // calculate the MP2 correlation
        if (hf.is_subcommand_used("mp2")) {
            // throw an error if no coulomb
            if (program.get<bool>("--no-coulomb")) throw std::runtime_error("I'm sorry, you need the coulomb tensor for the MP2 method.");

            // optimize the molecule with MP2 method
            if (mp2.is_used("-o")) {
                // print the MP2 optimization method header and optimize the molecule
                std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 OPTIMIZATION\n" << std::string(104, '-') << "\n\n";

                // perform the optimization
                data = MP(data).Optimizer.mp2();

                // print the optimized coordinates
                std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << data.system.coords << std::endl;

                // print the new distance matrix
                if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << data.system.dists << std::endl; 
            }

            // print the MP2 correlation method header
            std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 CORRELATION ENERGY\n" << std::string(104, '-') + "\n";

            // transform the coulomb tensor
            std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(data.intsmo.J = Transform::Coulomb(data.ints.J, data.roothaan.C))
            if (CONTAINS(mp2print, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << data.intsmo.J;} std::cout << "\n";
            if (CONTAINS(mp2save, "jmo") || CONTAINS(save, "all")) Eigen::Write("JMO.mat", data.intsmo.J);

            // do the calculation
            data = MP(data).mp2();

            // print the gradient and norm
            std::cout << "\nMP2 CORRELATION ENERGY: " << data.mp.Ecorr << std::endl << "FINAL ";
            std::cout << "MP2 ENERGY: " << data.roothaan.E + data.mp.Ecorr << std::endl;

            // calculate the MP2 nuclear gradient
            if (mp2.is_used("-g")) {
                // print the MP2 gradient method header and perform the calculation
                if (mp2.is_used("--numgrad")) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR MP2 METHOD\n" << std::string(104, '-');
                else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR MP2 METHOD\n" << std::string(104, '-') << "\n\n";

                // perform the calculation
                if (!mp2.is_used("-o")) data = MP(data).Gradient.mp2();

                // print the gradient results
                std::cout << data.mp.grad.G << "\n\nGRADIENT NORM: ";
                std::printf("%.2e\n", data.mp.grad.G.norm());
            }
        }

        // calculate ci correlation
        if (hf.is_subcommand_used("ci")) {
            // throw an error if no coulomb and print the CI method header
            if (program.get<bool>("--no-coulomb")) throw std::runtime_error("I'm sorry, you need the coulomb tensor for CI.");
            std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n" << std::string(104, '-') + "\n";

            // transform the coulomb tensor
            std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(data.intsmo.J = Transform::Coulomb(data.ints.J, data.roothaan.C))
            if (CONTAINS(ciprint, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << data.intsmo.J;} std::cout << "\n";
            if (CONTAINS(cisave, "jmo") || CONTAINS(save, "all")) Eigen::Write("JMO.mat", data.intsmo.J);

            // do the calculation
            data = CI(data).cid();

            // print and save the result matrices
            if (CONTAINS(ciprint, "cie") || CONTAINS(print, "all")) {std::cout << "\n" << data.ci.eig << "\n";}
            if (CONTAINS(ciprint, "cih") || CONTAINS(print, "all")) {std::cout << "\n" << data.ci.H << "\n";}
            if (CONTAINS(ciprint, "cic") || CONTAINS(print, "all")) {std::cout << "\n" << data.ci.C << "\n";}
            if (CONTAINS(cisave, "cie") || CONTAINS(save, "all")) Eigen::Write("CIE.mat", data.ci.eig);
            if (CONTAINS(cisave, "cih") || CONTAINS(save, "all")) Eigen::Write("CIH.mat", data.ci.H);
            if (CONTAINS(cisave, "cic") || CONTAINS(save, "all")) Eigen::Write("CIC.mat", data.ci.C);

            // print the gradient and norm
            std::cout << "\nCI CORRELATION ENERGY: " << data.ci.Ecorr << std::endl;
            std::cout << "FINAL CI ENERGY: " << data.roothaan.E + data.ci.Ecorr << std::endl;
        }
    }
}

Data Distributor::integrals(Data data) const {
    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(data.ints.S = Integral::Overlap(data.system))
    if (CONTAINS(print, "s") || CONTAINS(print, "all")) std::cout << "\n" << data.ints.S << std::endl;
    if (CONTAINS(save, "s") || CONTAINS(save, "all")) Eigen::Write("S.mat", data.ints.S);

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(data.ints.T = Integral::Kinetic(data.system))
    if (CONTAINS(print, "t") || CONTAINS(print, "all")) std::cout << "\n" << data.ints.T << std::endl;
    if (CONTAINS(save, "t") || CONTAINS(save, "all")) Eigen::Write("T.mat", data.ints.T);

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(data.ints.V = Integral::Nuclear(data.system))
    if (CONTAINS(print, "v") || CONTAINS(print, "all")) std::cout << "\n" << data.ints.V << std::endl;
    if (CONTAINS(save, "v") || CONTAINS(save, "all")) Eigen::Write("V.mat", data.ints.V);

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(data.ints.J = Integral::Coulomb(data.system))}
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(print, "j") || CONTAINS(print, "all"))) {std::cout << "\n" << data.ints.J;} std::cout << "\n";
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(save, "j") || CONTAINS(save, "all"))) Eigen::Write("J.mat", data.ints.J);

    // if derivatives of the integrals are needed
    if ((hf.is_used("-g") || hf.is_used("-o")) && !hf.is_used("numgrad")) {
        // calculate the overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(data.ints.dS = Integral::dOverlap(data.system))
        if (CONTAINS(print, "ds")) std::cout << "\n" << data.ints.dS << std::endl;

        // calculate the kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(data.ints.dT = Integral::dKinetic(data.system))
        if (CONTAINS(print, "dt")) std::cout << "\n" << data.ints.dT << std::endl;

        // calculate the nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(data.ints.dV = Integral::dNuclear(data.system))
        if (CONTAINS(print, "dv")) std::cout << "\n" << data.ints.dV << std::endl;

        // calculate the electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(data.ints.dJ = Integral::dCoulomb(data.system))
        if (CONTAINS(print, "dj")) {std::cout << "\n" << data.ints.dJ;} std::cout << "\n";
    }

    // return integrals
    return data;
}
