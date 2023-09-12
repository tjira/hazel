#include "distributor.h"

#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // initialize general subcommands
    ints = argparse::ArgumentParser("ints", "0.1", argparse::default_arguments::none), hf = argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none);
    mp2 = argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none), ci = argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none);
    qd = argparse::ArgumentParser("qd", "0.1", argparse::default_arguments::none), md = argparse::ArgumentParser("md", "0.1", argparse::default_arguments::none);

    // initialize MD subcommands
    mdhf = argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none), mdmp2 = argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none);

    // add positional arguments to the main argument parser
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-c", "--charge").help("-- Charge of the system.").default_value(0).scan<'i', int>();
    program.add_argument("-f", "--file").help("-- Path to the system to use in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads to use.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("-s", "--spin").help("-- Spin multiplicity of the system.").default_value(1).scan<'i', int>();
    program.add_argument("--no-center").help("-- Disable the molecule centering.").default_value(false).implicit_value(true);
    program.add_argument("--no-coulomb").help("-- Disable calculation of the coulomb tensor.").default_value(false).implicit_value(true);
    program.add_argument("--show-bases").help("-- Print all the available bases.").default_value(false).implicit_value(true);

    // add positional arguments to the integral argument parser
    ints.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    ints.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the HF argument parser
    hf.add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    hf.add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    hf.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    hf.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    hf.add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    hf.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    hf.add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    mp2.add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    mp2.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    mp2.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    mp2.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    mp2.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the CI argument parser
    ci.add_argument("-e", "--excitations").help("-- Excitations to consider.").default_value("d");
    ci.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    ci.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the QD argument parser
    qd.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    qd.add_argument("-s", "--step").help("-- Dynamics time step in atomic units.").default_value(0.1).scan<'g', double>();
    qd.add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(1000).scan<'i', int>();
    qd.add_argument("-r", "--range").help("-- Grid range in all dimensions.").default_value(16.0).scan<'g', double>();
    qd.add_argument("-n", "--nstate").help("-- Number of states to consider.").default_value(3).scan<'i', int>();
    qd.add_argument("-p", "--points").help("-- Number of points on the grid.").default_value(1024).scan<'i', int>();
    qd.add_argument("-t", "--thresh").help("-- Threshold for conververgence in ITP loop.").default_value(1e-8).scan<'g', double>();
    qd.add_argument("--no-real").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the MD argument parser
    md.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    md.add_argument("-s", "--step").help("-- Dynamics time step in atomic units.").default_value(0.5).scan<'g', double>();
    md.add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(100).scan<'i', int>();
    md.add_argument("-o", "--output").help("-- Output of the trajectory.").default_value("trajectory.xyz");

    // add positional arguments to the MD HF argument parser
    mdhf.add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    mdhf.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    mdhf.add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    mdhf.add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    mdhf.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD MP2 argument parser
    mdmp2.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    mdmp2.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add the parsers
    program.add_subparser(ints), program.add_subparser(hf), hf.add_subparser(mp2), hf.add_subparser(ci);
    program.add_subparser(qd), program.add_subparser(md), md.add_subparser(mdhf), mdhf.add_subparser(mdmp2);

    // parse the arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help if requested
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("ints") && ints.get<bool>("-h")) {
        std::cout << ints.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && hf.get<bool>("-h")) {
        std::cout << hf.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && hf.is_subcommand_used("mp2") && mp2.get<bool>("-h")) {
        std::cout << mp2.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && hf.is_subcommand_used("ci") && ci.get<bool>("-h")) {
        std::cout << ci.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("qd") && qd.get<bool>("-h")) {
        std::cout << qd.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && md.get<bool>("-h")) {
        std::cout << md.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && md.is_subcommand_used("hf") && mdhf.get<bool>("-h")) {
        std::cout << mdhf.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && md.is_subcommand_used("hf") && mdhf.is_subcommand_used("mp2") && mdmp2.get<bool>("-h")) {
        std::cout << mdmp2.help().str(); exit(EXIT_SUCCESS);
    }

    // print all tha available bases if requested
    if (program.get<bool>("--show-bases")) {
        for (const auto& file : std::filesystem::directory_iterator(std::string(DATADIR) + "/basis")) {
            if (file.path().filename() != "basis.sh") std::cout << file.path().filename().replace_extension().c_str() << std::endl;
        }
        exit(EXIT_SUCCESS);
    }

    // set the path to the basis functions
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        std::string path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
    }
    
    // set number of threads to use and cout flags
    std::cout << std::fixed << std::setprecision(14);
    nthread = program.get<int>("--nthread");
}

Distributor::~Distributor() {
    std::cout << "\n" + std::string(104, '-') << std::endl;
    std::printf("TOTAL EXECUTION TIME: %s\n", Timer::Format(Timer::Elapsed(start)).c_str());
    std::cout << std::string(104, '-') << std::endl;
}

void Distributor::run() {
    // initialize the system and extract general printing flag
    system = System(program.get("-f"), program.get("-b"), program.get<int>("-c"), program.get<int>("-s"));
    print = program.get<std::vector<std::string>>("-p");

    // print the title with number of threads
    std::cout << "QUANTUM HAZEL" << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, Timer::Local().c_str());
    std::printf("\nCOMPILER VERSION AND FLAGS: MACHINE %s, GCC %d.%d.%d (%s)", []{utsname unr; uname(&unr); return unr;}().machine, __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d, LIBXC %s\n", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION, XC_VERSION);
    std::printf("\nAVAILABLE CORES: %d\nUSED THREADS: %d\n", std::thread::hardware_concurrency(), nthread);

    // print the system block
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + system.basis + ")\n" << std::string(104, '-') + "\n\n";
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);
    std::cout << "\nSYSTEM COORDINATES\n" << system.coords << std::endl; 

    // center the molecule if not disabled
    if (!program.get<bool>("--no-center")) {
        Matrix dir(system.atoms.size(), 3); dir.rowwise() -= system.coords.colwise().sum() / system.atoms.size();
        system.move(dir * A2BOHR); std::cout << "\nCENTERED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    }

    // print the distances if requested
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nDISTANCE MATRIX\n" << system.dists << std::endl; 

    // extract printing options
    if (program.is_subcommand_used("ints")) {
        intsprint = ints.get<std::vector<std::string>>("-p");
    } else if (program.is_subcommand_used("hf")) {
        hfprint = hf.get<std::vector<std::string>>("-p");
        if (hf.is_subcommand_used("mp2")) {
            mp2print = mp2.get<std::vector<std::string>>("-p");
        } else if (hf.is_subcommand_used("ci")) {
            ciprint = ci.get<std::vector<std::string>>("-p");
        }
    }

    // extract HF options
    if (program.is_subcommand_used("hf")) {
        if (program.get<int>("-s") == 1) {
            rhfopt.diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
            rhfopt.maxiter = hf.get<int>("-i"), rhfopt.thresh = hf.get<double>("-t");
            rhfopt.nocoulomb = program.get<bool>("--no-coulomb");
        } else {
            uhfopt.diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
            uhfopt.maxiter = hf.get<int>("-i"), uhfopt.thresh = hf.get<double>("-t");
            uhfopt.nocoulomb = program.get<bool>("--no-coulomb");
        }
    } else if (program.is_subcommand_used("md") && md.is_subcommand_used("hf")) {
        if (program.get<int>("-s") == 1) {
            rhfopt.diis = {mdhf.get<std::vector<int>>("-d").at(0), mdhf.get<std::vector<int>>("-d").at(1)};
            rhfopt.maxiter = mdhf.get<int>("-i"), rhfopt.thresh = mdhf.get<double>("-t");
            rhfopt.nocoulomb = program.get<bool>("--no-coulomb");
        } else {
            uhfopt.diis = {mdhf.get<std::vector<int>>("-d").at(0), mdhf.get<std::vector<int>>("-d").at(1)};
            uhfopt.maxiter = mdhf.get<int>("-i"), uhfopt.thresh = mdhf.get<double>("-t");
            uhfopt.nocoulomb = program.get<bool>("--no-coulomb");
        }
    }

    // transform the print vectors to lowercase
    for (auto& element : mp2print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : hfprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : ciprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // calculate the integrals if needed
    if (program.is_subcommand_used("ints") || program.is_subcommand_used("hf")) integrals();

    // optimize the molecule
    if (program.is_subcommand_used("hf")) {
        if (hf.is_used("-o")) {
            if (program.get<int>("-s") == 1) rhfo();
            else throw std::runtime_error("OPTIMIZATION FOR UHF NOT IMPLEMENTED");
        }
        else if (hf.is_subcommand_used("mp2")) {
            if (mp2.is_used("-o")) {
                if (program.get<int>("-s") == 1) rmp2o();
                else throw std::runtime_error("OPTIMIZATION FOR UMP2 NOT IMPLEMENTED");
            }
        }
    }

    // distribute the calculations
    if (program.is_subcommand_used("md")) dynamics();
    if (program.is_subcommand_used("qd")) qdyn();
    if (program.is_subcommand_used("hf")) {
        if (program.get<int>("-s") == 1) rhfrun();
        else uhfrun();
    }
}

void Distributor::rcirun() {
    // print the CI method header and define J in MO basis
    std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo; CI::ResultsRestricted rcires;

    // transform the coulomb tensor
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (CONTAINS(ciprint, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << Jmo;} std::cout << "\n";

    // do the calculation
    if (ci.get("-e") == "s") rcires = CI({rhfres}).rcis(system, Jmo);
    if (ci.get("-e") == "d") rcires = CI({rhfres}).rcid(system, Jmo);

    // print the result matrices
    if (CONTAINS(ciprint, "cih") || CONTAINS(print, "all")) std::cout << "\nCI HAMILTONIAN\n" << rcires.H << "\n";
    if (CONTAINS(ciprint, "cie") || CONTAINS(print, "all")) std::cout << "\nCI ENERGIES\n" << Matrix(rcires.eig) << "\n";
    if (CONTAINS(ciprint, "cic") || CONTAINS(print, "all")) std::cout << "\nCI EXPANSION COEFFICIENTS\n" << rcires.C << "\n";

    // print the gradient and norm
    std::cout << "\nCI CORRELATION ENERGY: " << rcires.Ecorr << std::endl;
    std::cout << "FINAL CI ENERGY: " << rhfres.E + rcires.Ecorr << std::endl;
}

void Distributor::rhfrun() {
    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", rhfopt.maxiter, rhfopt.thresh, rhfopt.diis.start, rhfopt.diis.keep);

    // perform the RHF calculation
    rhfres = HF(rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print the resulting RHF matrices and energies
    if (CONTAINS(hfprint, "eps") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(rhfres.eps) << std::endl;
    if (CONTAINS(hfprint, "c") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX\n" << rhfres.C << std::endl;
    if (CONTAINS(hfprint, "d") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX\n" << rhfres.D << std::endl;
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << rhfres.E << std::endl;

    // gradient and frequency
    if (hf.is_used("-g")) rhfg();
    if (hf.is_used("-f")) rhff();

    // post RHF methods
    if (hf.is_subcommand_used("mp2")) rmp2run();
    if (hf.is_subcommand_used("ci")) rcirun();
}

void Distributor::rhff() {
    // print the Hessian header
    if (hf.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H;

    // perform the Hessian calculation
    if (hf.get<std::vector<double>>("-f").at(0)) H = Hessian({hf.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RHF NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({hf.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS\n" << std::string(104, '-');
    std::cout << "\n\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::rhfg() {
    // print the header for gradient calculation and define the gradient matrix
    if (hf.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix G; 

    // calculate the numerical or analytical gradient
    if (hf.get<std::vector<double>>("-g").at(0)) G = Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else G = Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, rhfres);

    // print the gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rhfo() {
    // print the RHF optimization header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n\n";

    // perform the optimization
    system = Optimizer({hf.get<double>("-o")}).optimize(system, Lambda::EGHF(rhfopt, hf.get<std::vector<double>>("-g"), rhfres.D));

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::uhfrun() {
    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nUNRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", uhfopt.maxiter, uhfopt.thresh, uhfopt.diis.start, uhfopt.diis.keep);

    // perform the RHF calculation
    uhfres = HF(uhfopt).uscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print the resulting RHF matrices and energies
    if (CONTAINS(hfprint, "epsa") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES (ALPHA)\n" << Matrix(uhfres.epsa) << std::endl;
    if (CONTAINS(hfprint, "epsb") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES (BETA)\n" << Matrix(uhfres.epsb) << std::endl;
    if (CONTAINS(hfprint, "ca") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX (ALPHA)\n" << uhfres.Ca << std::endl;
    if (CONTAINS(hfprint, "cb") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX (BETA)\n" << uhfres.Cb << std::endl;
    if (CONTAINS(hfprint, "da") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX (ALPHA)\n" << uhfres.Da << std::endl;
    if (CONTAINS(hfprint, "db") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX (BETA)\n" << uhfres.Db << std::endl;
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << uhfres.E << std::endl;

    // gradient and frequency
    if (hf.is_used("-g")) uhfg();
    if (hf.is_used("-f")) uhff();

    // post UHF methods
    if (hf.is_subcommand_used("mp2")) throw std::runtime_error("UMP2 NOT IMPLEMENTED");
    if (hf.is_subcommand_used("ci")) throw std::runtime_error("UCI NOT IMPLEMENTED");
}

void Distributor::uhff() {
    // print the Hessian header
    if (hf.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H;

    // perform the Hessian calculation
    if (hf.get<std::vector<double>>("-f").at(0)) H = Hessian({hf.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR UHF NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({hf.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS\n" << std::string(104, '-');
    std::cout << "\n\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::uhfg() {
    // print the header for gradient calculation and define the gradient matrix
    if (hf.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix G; 

    // calculate the numerical or analytical gradient
    if (hf.get<std::vector<double>>("-g").at(0)) G = Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");

    // print the gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rmp2run() {
    // print the MP2 method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;

    // transform the coulomb tensor to MO basis
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (CONTAINS(mp2print, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << Jmo;} std::cout << "\n";

    // perform the MP2 calculation
    double Ecorr = MP({rhfres}).rmp2(system, Jmo);

    // print the gradient and it's norm
    std::cout << "\nMP2 CORRELATION ENERGY: " << Ecorr << std::endl << "FINAL ";
    std::cout << "MP2 ENERGY: " << rhfres.E + Ecorr << std::endl;

    // gradient and frequency
    if (mp2.is_used("-g")) rmp2g();
    if (mp2.is_used("-f")) rmp2f();
}

void Distributor::rmp2f() {
    // print the MP2 frequency calculation header
    if (mp2.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H; 

    // perform the Hessian calculation
    if (mp2.get<std::vector<double>>("-f").at(0)) H = Hessian({mp2.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RMP2 NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({mp2.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the MP2 frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 FREQUENCY ANALYSIS\n" << std::string(104, '-') + "\n";
    std::cout << "\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::rmp2g() {
    // print the MP2 gradient method header
    if (mp2.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n"; Matrix G;

    // perform the MP2 gradient calculation
    if (mp2.get<std::vector<double>>("-g").at(0)) G = Gradient({mp2.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RMP2 NOT IMPLEMENTED");

    // print the MP2 gradient with its norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rmp2o() {
    // print the MP2 optimization method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 OPTIMIZATION\n" << std::string(104, '-') << "\n\n";

    // perform the optimization
    system = Optimizer({mp2.get<double>("-o")}).optimize(system, Lambda::EGMP2(rhfopt, mp2.get<std::vector<double>>("-g"), rhfres.D));

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl;
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::qdyn() {
    // print the dynamics header
    std::cout << "\n" + std::string(104, '-') + "\nQUANTUM DYNAMICS\n" << std::string(104, '-') + "\n\n";

    // perform the dynamics
    Qdyn::Results qdres = Qdyn({qd.get<int>("-p"), qd.get<int>("-i"), qd.get<int>("-n"), qd.get<double>("-r"), qd.get<double>("-s"), qd.get<double>("-t"), qd.get<bool>("--no-real")}).run(system);

    // print the energies
    if (qd.get<bool>("--no-real")) std::cout << "\nIMAGINARY TIME PROPAGATION ENERGIES\n" << Matrix(qdres.energy) << std::endl;

    // save the wavefunction
    Qdyn::wfnsave(qdres.r, qdres.states, "wavefunction.dat");
}

void Distributor::dynamics() {
    // print the dynamics header
    std::cout << "\n" + std::string(104, '-') + "\nMOLECULAR DYNAMICS\n" << std::string(104, '-') + "\n\n";

    // define the anonymous function for gradient
    std::function<std::tuple<double, Matrix>(System&)> egfunc;

    if (md.is_subcommand_used("hf") && program.get<int>("-s") != 1) {
        throw std::runtime_error("DYNAMICS NOT IMPLEMENTED FOR UNRESTRICTED CALCULATION");
    }

    // get the enrgy and gradient function
    if (md.is_subcommand_used("hf")) {
        egfunc = Lambda::EGHF(rhfopt, mdhf.get<std::vector<double>>("-g"), rhfres.D);
        if (mdhf.is_subcommand_used("mp2")) {
            egfunc = Lambda::EGMP2(rhfopt, mdmp2.get<std::vector<double>>("-g"), rhfres.D);
        }
    } else throw std::runtime_error("INVALID METHOD FOR DYNAMICS SPECIFIED");

    // perform the dynamics
    Dynamics({md.get<int>("-i"), md.get<double>("-s"), md.get("-o")}).run(system, egfunc);
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.S);
    if (CONTAINS(intsprint, "s") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.ints.S << std::endl;
    if (CONTAINS(hfprint, "s") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.ints.S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.T);
    if (CONTAINS(intsprint, "t") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.ints.T << std::endl;
    if (CONTAINS(hfprint, "t") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.ints.T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.V);
    if (CONTAINS(intsprint, "v") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.ints.V << std::endl;
    if (CONTAINS(hfprint, "v") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.ints.V << std::endl;

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    if (!program.get<bool>("--no-coulomb")) std::cout << " " << Eigen::MemTensor(system.ints.J);
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(intsprint, "j") || CONTAINS(intsprint, "all"))) {std::cout << "\n" << system.ints.J;}
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(hfprint, "j") || CONTAINS(hfprint, "all"))) {std::cout << "\n" << system.ints.J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.is_subcommand_used("ints") || ((hf.is_used("-g") || hf.is_used("-o")) && !hf.get<std::vector<double>>("-g").at(0))) {
        // calculate the derivative of overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dS);
        if (CONTAINS(intsprint, "ds") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.dints.dS << std::endl;
        if (CONTAINS(hfprint, "ds") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the derivative of kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dT);
        if (CONTAINS(intsprint, "dt") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.dints.dT << std::endl;
        if (CONTAINS(hfprint, "dt") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the derivative of nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dV);
        if (CONTAINS(intsprint, "dv") || CONTAINS(intsprint, "all")) std::cout << "\n" << system.dints.dV << std::endl;
        if (CONTAINS(hfprint, "dv") || CONTAINS(hfprint, "all")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the derivative of electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dJ);
        if (CONTAINS(intsprint, "dj") || CONTAINS(intsprint, "all")) {std::cout << "\n" << system.dints.dJ;}
        if (CONTAINS(hfprint, "dj") || CONTAINS(hfprint, "all")) {std::cout << "\n" << system.dints.dJ;} std::cout << "\n";
    }
}
