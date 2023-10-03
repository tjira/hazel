#include "distributor.h"

#define CONTAINS(V, E) [](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E)
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}
#define TOUPPER(S) [](std::string str) {for (auto& c : str) c = toupper(c); return str;}(S)

Distributor::Distributor(int argc, char** argv) : parsers(11), program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // add level 1 parsers
    parsers.push_back(argparse::ArgumentParser("ints", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("scan", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("qd", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("md", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));

    // add HF parsers
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));

    // add level 2 parsers
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").add_subparser(parsers.at(parsers.size() - 1));

    // add level 3 parsers
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    
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
    program.at<argparse::ArgumentParser>("ints").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("ints").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the HF argument parser
    program.at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-8).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the CI argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-e", "--excitations").help("-- Excitations to consider.").default_value("d");
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the SCAN argument parser
    program.at<argparse::ArgumentParser>("scan").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").add_argument("-o", "--output").help("-- Output of the PES energies.").default_value("pes.dat");

    // add positional arguments to the SCAN HF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add positional arguments to the SCAN MP2 argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the QD argument parser
    program.at<argparse::ArgumentParser>("qd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("qd").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.1).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(1000).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-n", "--nstate").help("-- Number of states to consider.").default_value(3).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-o", "--output").help("-- Output of the wavefunction.").default_value("wavefunction.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-f", "--potfile").help("-- File with the PES.").default_value("pes.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-t", "--thresh").help("-- Threshold for conververgence in ITP loop.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("--no-real").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the MD argument parser
    program.at<argparse::ArgumentParser>("md").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").add_argument("-o", "--output").help("-- Output of the trajectory.").default_value("trajectory.xyz");

    // add positional arguments to the MD HF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD MP2 argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // parse the arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help if requested
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("ints") && program.at<argparse::ArgumentParser>("ints").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("ints").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("hf").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("hf").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2") && program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci") && program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("scan") && program.at<argparse::ArgumentParser>("scan").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("scan").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("scan") && program.at<argparse::ArgumentParser>("scan").is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("scan") && program.at<argparse::ArgumentParser>("scan").is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2") && program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("qd") && program.at<argparse::ArgumentParser>("qd").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("qd").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && program.at<argparse::ArgumentParser>("md").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("md").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && program.at<argparse::ArgumentParser>("md").is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").help().str(); exit(EXIT_SUCCESS);
    } else if (program.is_subcommand_used("md") && program.at<argparse::ArgumentParser>("md").is_subcommand_used("hf") && program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2") && program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").get<bool>("-h")) {
        std::cout << program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").help().str(); exit(EXIT_SUCCESS);
    }

    // print all tha available bases if requested
    if (program.get<bool>("--show-bases")) {
        for (const auto& file : std::filesystem::directory_iterator(std::string(DATADIR) + "/basis")) {
            if (file.path().filename() != "basis.sh") std::cout << file.path().filename().replace_extension() << std::endl;
        }
        exit(EXIT_SUCCESS);
    }

    // set the path to the basis functions
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        auto path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path(); char buffer[4096];
        #ifdef _WIN32
        wcstombs(buffer, path.c_str(), path.string().size()); _putenv_s("LIBINT_DATA_PATH", buffer);
        #else
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
        #endif
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
    // initialize the system
    std::string basis = program.get("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
    std::ifstream stream(program.get("-f")); system = System(stream, basis, program.get<int>("-c"), program.get<int>("-s")); stream.close();

    // print the title with number of threads
    std::cout << "QUANTUM HAZEL" << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, Timer::Local().c_str());
    std::printf("\nCOMPILER VERSION AND FLAGS: MACHINE ---, GCC %d.%d.%d (%s)", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d\n", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
    std::printf("\nAVAILABLE CORES: %d\nUSED THREADS: %d\n", std::thread::hardware_concurrency(), nthread);

    // print the system block
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + program.get("-b") + ")\n" << std::string(104, '-') + "\n\n";
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);
    std::cout << "\nSYSTEM COORDINATES\n" << system.coords << std::endl; 

    // center the molecule if not disabled
    if (!program.get<bool>("--no-center")) {
        Matrix dir(system.atoms.size(), 3); dir.rowwise() -= system.coords.colwise().sum() / system.atoms.size();
        system.move(dir * A2BOHR); std::cout << "\nCENTERED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    }

    // print the distances if requested
    if (CONTAINS(program.get<std::vector<std::string>>("-p"), "dist")) std::cout << "\nDISTANCE MATRIX\n" << system.dists << std::endl; 

    // calculate the integrals if needed
    if (program.is_subcommand_used("ints") || program.is_subcommand_used("hf")) integrals();

    // optimize the molecule
    if (program.is_subcommand_used("hf")) {
        if (program.at<argparse::ArgumentParser>("hf").is_used("-o")) {
            if (program.get<int>("-s") == 1) rhfo(program.at<argparse::ArgumentParser>("hf"));
            else uhfo(program.at<argparse::ArgumentParser>("hf"));
        }
        else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").is_used("-o")) {
                if (program.get<int>("-s") == 1) rmp2o(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2"));
                else throw std::runtime_error("OPTIMIZATION FOR UMP2 NOT IMPLEMENTED");
            }
        }
    }

    // distribute the calculations
    if (program.is_subcommand_used("scan")) scan(program.at<argparse::ArgumentParser>("scan"));
    if (program.is_subcommand_used("md")) dynamics(program.at<argparse::ArgumentParser>("md"));
    if (program.is_subcommand_used("qd")) qdyn(program.at<argparse::ArgumentParser>("qd"));
    if (program.is_subcommand_used("hf")) {
        if (program.get<int>("-s") == 1) rhfrun(program.at<argparse::ArgumentParser>("hf"));
        else uhfrun(program.at<argparse::ArgumentParser>("hf"));
    }
}

void Distributor::rcirun(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the CI method header and define J in MO basis
    std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;
    std::printf("\n-- EXCITATIONS: %s\n", TOUPPER(parser.get("-e")).c_str()); CI::ResultsRestricted rcires;

    // transform the coulomb tensor
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "jmo")) {std::cout << "\n" << Jmo;} std::cout << "\n";

    // do the calculation
    if (parser.get("-e") == "s") rcires = CI({rhfres}).rcis(system, Jmo);
    if (parser.get("-e") == "d") rcires = CI({rhfres}).rcid(system, Jmo);

    // print the result matrices
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "cih")) std::cout << "\nCI HAMILTONIAN\n" << rcires.H << "\n";
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "cie")) std::cout << "\nCI ENERGIES\n" << Matrix(rcires.eig) << "\n";
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "cic")) std::cout << "\nCI EXPANSION COEFFICIENTS\n" << rcires.C << "\n";

    // print the gradient and norm
    std::cout << "\nCI CORRELATION ENERGY: " << rcires.Ecorr << std::endl;
    std::cout << "FINAL CI ENERGY: " << rhfres.E + rcires.Ecorr << std::endl;
}

void Distributor::rhfrun(argparse::ArgumentParser& parser) {
    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", rhfopt.maxiter, rhfopt.thresh, rhfopt.diis.start, rhfopt.diis.keep);

    // perform the RHF calculation
    auto rhfres = HF(rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print the resulting RHF matrices and energies
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "eps")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(rhfres.eps) << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "c")) std::cout << "\nCOEFFICIENT MATRIX\n" << rhfres.C << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "d")) std::cout << "\nDENSITY MATRIX\n" << rhfres.D << std::endl;
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << rhfres.E << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) rhfg(parser, rhfres);
    if (parser.is_used("-f")) rhff(parser, rhfres);

    // post RHF methods
    if (parser.is_subcommand_used("mp2")) rmp2run(parser.at<argparse::ArgumentParser>("mp2"), rhfres);
    if (parser.is_subcommand_used("ci")) rcirun(parser.at<argparse::ArgumentParser>("ci"), rhfres);
}

void Distributor::rhff(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the Hessian header
    if (parser.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H;

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RHF NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS\n" << std::string(104, '-');
    std::cout << "\n\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::rhfg(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix G; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, rhfres);

    // print the gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rhfo(argparse::ArgumentParser& parser) {
    // print the RHF optimization header and extract the HF options
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n\n";
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGHF(rhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "dist")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::uhfrun(argparse::ArgumentParser& parser) {
    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nUNRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", uhfopt.maxiter, uhfopt.thresh, uhfopt.diis.start, uhfopt.diis.keep);

    // perform the RHF calculation
    auto uhfres = HF(uhfopt).uscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print the resulting RHF matrices and energies
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "epsa")) std::cout << "\nORBITAL ENERGIES (ALPHA)\n" << Matrix(uhfres.epsa) << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "epsb")) std::cout << "\nORBITAL ENERGIES (BETA)\n" << Matrix(uhfres.epsb) << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "ca")) std::cout << "\nCOEFFICIENT MATRIX (ALPHA)\n" << uhfres.Ca << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "cb")) std::cout << "\nCOEFFICIENT MATRIX (BETA)\n" << uhfres.Cb << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "da")) std::cout << "\nDENSITY MATRIX (ALPHA)\n" << uhfres.Da << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "db")) std::cout << "\nDENSITY MATRIX (BETA)\n" << uhfres.Db << std::endl;
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << uhfres.E << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) uhfg(parser, uhfres);
    if (parser.is_used("-f")) uhff(parser, uhfres);

    // post UHF methods
    if (parser.is_subcommand_used("mp2")) throw std::runtime_error("UMP2 NOT IMPLEMENTED");
    if (parser.is_subcommand_used("ci")) throw std::runtime_error("UCI NOT IMPLEMENTED");
}

void Distributor::uhff(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres) {
    // print the Hessian header
    if (parser.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H;

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR UHF NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS\n" << std::string(104, '-');
    std::cout << "\n\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::uhfg(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix G; 

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");

    // print the gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::uhfo(argparse::ArgumentParser& parser) {
    // print the UHF optimization header and extract the HF options
    std::cout << "\n" + std::string(104, '-') + "\nUNRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n\n";
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGHF(uhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "dist")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::rmp2run(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the MP2 method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;

    // transform the coulomb tensor to MO basis
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "jmo")) {std::cout << "\n" << Jmo;} std::cout << "\n";

    // perform the MP2 calculation
    double Ecorr = MP({rhfres}).rmp2(system, Jmo);

    // print the gradient and it's norm
    std::cout << "\nMP2 CORRELATION ENERGY: " << Ecorr << std::endl << "FINAL ";
    std::cout << "MP2 ENERGY: " << rhfres.E + Ecorr << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) rmp2g(parser, rhfres);
    if (parser.is_used("-f")) rmp2f(parser, rhfres);
}

void Distributor::rmp2f(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the MP2 frequency calculation header
    if (parser.get<std::vector<double>>("-f").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    std::cout << std::string(104, '-') + "\n\n"; Matrix H; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RMP2 NOT IMPLEMENTED");

    // print the Hessian with its norm and perform the frequency calculation
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the MP2 frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 FREQUENCY ANALYSIS\n" << std::string(104, '-') + "\n";
    std::cout << "\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::rmp2g(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the MP2 gradient method header and extract the HF options
    if (parser.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";
    else {std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";} Matrix G;
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the MP2 gradient calculation
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RMP2 NOT IMPLEMENTED");

    // print the MP2 gradient with its norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rmp2o(argparse::ArgumentParser& parser) {
    // extract the HF options and print the MP2 optimization method header
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 OPTIMIZATION\n" << std::string(104, '-') << "\n\n";

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGMP2(rhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl;
    if (CONTAINS(parser.get<std::vector<std::string>>("-p"), "dist")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::scan(argparse::ArgumentParser& parser) {
    // print the energy scan header
    std::cout << "\n" + std::string(104, '-') + "\nMOVIE ENERGY SCAN\n" << std::string(104, '-') + "\n";
    std::printf("\nITER       Eel [Eh]           TIME    \n");

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // define the anonymous function for gradient
    std::function<double(System)> efunc; std::vector<double> energies;

    // get the energy and gradient function
    if (parser.is_subcommand_used("hf") && program.get<int>("-s") == 1) {
        efunc = Lambda::EHF(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            efunc = Lambda::EMP2(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.is_subcommand_used("hf") && program.get<int>("-s") != 1) {
        efunc = Lambda::EHF(uhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            throw std::runtime_error("UMP2 NOT IMPLEMENTED");
        }
    } else throw std::runtime_error("INVALID METHOD FOR SCAN SPECIFIED");

    // count the number of geometries
    std::ifstream stream(program.get("-f")); std::string line; int lines;
    for (lines = 0; std::getline(stream, line); lines++);
    int geoms = lines / (system.atoms.size() + 2);
    stream.clear(), stream.seekg(0);

    // perform the scan
    for (int i = 0; i < geoms; i++) {
        // create the system
        std::string basis = program.get("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
        system = System(stream, basis, program.get<int>("-c"), program.get<int>("-s")); Timer::Timepoint start = Timer::Now();

        // start the timer and calculate the energy
        double E = efunc(system); energies.push_back(E);

        // print the energy
        std::printf("%4d %20.14f %s\n", i + 1, E, Timer::Format(Timer::Elapsed(start)).c_str());
    }
    
    // save the file
    std::ofstream file(parser.get("-o"));
    file << std::fixed << std::setprecision(14) << "# i              E\n";
    for (size_t i = 0; i < energies.size(); i++) {
        file << std::setw(5) << i << " " << std::setw(20) << energies.at(i) << "\n";
    }
}

void Distributor::dynamics(argparse::ArgumentParser& parser) {
    // print the dynamics header
    std::cout << "\n" + std::string(104, '-') + "\nMOLECULAR DYNAMICS\n" << std::string(104, '-') + "\n\n";

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // define the anonymous function for gradient
    std::function<std::tuple<double, Matrix>(System&)> egfunc;

    // get the energy and gradient function
    if (parser.is_subcommand_used("hf") && program.get<int>("-s") == 1) {
        egfunc = Lambda::EGHF(rhfopt, parser.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            egfunc = Lambda::EGMP2(rhfopt, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.is_subcommand_used("hf") && program.get<int>("-s") != 1) {
        egfunc = Lambda::EGHF(uhfopt, parser.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            throw std::runtime_error("UMP2 NOT IMPLEMENTED");
        }
    } else throw std::runtime_error("INVALID METHOD FOR DYNAMICS SPECIFIED");

    // perform the dynamics
    MD({parser.get<int>("-i"), parser.get<double>("-s"), parser.get("-o")}).run(system, egfunc);
}

void Distributor::qdyn(argparse::ArgumentParser& parser) {
    // print the dynamics header
    std::cout << "\n" + std::string(104, '-') + "\nQUANTUM DYNAMICS\n" << std::string(104, '-') + "\n\n";

    // perform the dynamics
    QD::Results qdres = QD({parser.get("-f"), parser.get<int>("-i"), parser.get<int>("-n"), parser.get<double>("-s"), parser.get<double>("-t"), parser.get<bool>("--no-real")}).run(system);

    // print the energies
    if (parser.get<bool>("--no-real")) std::cout << "\nIMAGINARY TIME PROPAGATION ENERGIES\n" << Matrix(qdres.energy) << std::endl;

    // save the wavefunction
    Utility::SaveWavefunction(parser.get<std::string>("-o"), qdres.r, qdres.states, qdres.energy);
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.S);
    if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;
    if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.T);
    if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;
    if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.V);
    if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;
    if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    if (!program.get<bool>("--no-coulomb")) std::cout << " " << Eigen::MemTensor(system.ints.J);
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;

    // print new line
    std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.is_subcommand_used("ints") || ((program.at<argparse::ArgumentParser>("hf").is_used("-g") || program.at<argparse::ArgumentParser>("hf").is_used("-o")) && !program.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g").at(0))) {
        // calculate the derivative of overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dS);
        if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;
        if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the derivative of kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dT);
        if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;
        if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the derivative of nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dV);
        if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;
        if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the derivative of electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dJ);
        if (program.is_subcommand_used("ints") && CONTAINS(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;
        if (program.is_subcommand_used("hf") && CONTAINS(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;

        // print new line
        std::cout << "\n";
    }
}
