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
    program.add_argument("-c", "--charge").help("-- Molecular charge.").default_value(0).scan<'i', int>();
    program.add_argument("-f", "--file").help("-- Quantum system to use in .xyz file format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads to use.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("--center").help("-- Center the molecule before doing any calculation.").default_value(false).implicit_value(true);
    program.add_argument("--no-coulomb").help("-- Disable calculation of the coulomb tensor.").default_value(false).implicit_value(true);

    // add positional arguments to the HF argument parser
    hf.add_argument("-d", "--diis").help("-- Start iteration and Fock history length for DIIS.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    hf.add_argument("-f", "--frequency").help("-- Enable analytical (0) or numerical (1) frequency calculation.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    hf.add_argument("-g", "--gradient").help("-- Enable analytical (0) or numerical (1) gradient calculation.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    hf.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    hf.add_argument("-m", "--maxiter").help("-- Maximum number of iterations to do in iterative calculations.").default_value(100).scan<'i', int>();
    hf.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    hf.add_argument("-t", "--thresh").help("-- Threshold for conververgence.").default_value(1e-12).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    mp2.add_argument("-f", "--frequency").help("-- Enable analytical (0) or numerical (1) frequency calculation.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    mp2.add_argument("-g", "--gradient").help("-- Enable analytical (0) or numerical (1) gradient calculation.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    mp2.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    mp2.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    mp2.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments tinput001.outo the CI argument parser
    ci.add_argument("-e", "--excitations").help("-- Define what excitations to use in the CI calculation.").default_value("s");
    ci.add_argument("-h", "--help").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    ci.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

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
    // initialize the system and extract flags
    system = System(program.get("-f"), program.get("-b"), program.get<int>("-c"), 1);
    print = program.get<std::vector<std::string>>("-p");

    // check if unrestricted calculation needed
    if (system.charge % 2) throw std::runtime_error("SPIN UNRESTRICTED CALCULATIONS ARE NOT SUPPORTED YET");

    // print the title with number of threads
    std::cout << "QUANTUM HAZEL" << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, LOCALTIME.c_str());
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
    std::printf("\nCOMPILER VERSION: GCC %d.%d.%d\nCOMPILER FLAGS: %s\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nAVAILABLE CORES: %d\nUSED THREADS: %d\n", std::thread::hardware_concurrency(), nthread);

    // print the system block
    std::cout << "\n" + std::string(104, '-') + "\nSYSTEM SPECIFICATION (" + system.basis + ")\n" << std::string(104, '-') + "\n\n";
    std::printf("-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
    std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);
    std::cout << "\nSYSTEM COORDINATES\n" << system.coords << std::endl; 

    // center the molecule if requested
    if (program.get<bool>("--center")) {
        Matrix dir(system.atoms.size(), 3); dir.rowwise() -= system.coords.colwise().sum() / system.atoms.size();
        system.move(dir * A2BOHR); std::cout << "\nCENTERED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    }

    // print the distances if requested
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nDISTANCE MATRIX\n" << system.dists << std::endl; 

    // extract all the options
    if (program.is_subcommand_used("hf")) {
        // printing options
        hfprint = hf.get<std::vector<std::string>>("-p");

        // scf options
        rhfopt.diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
        rhfopt.maxiter = hf.get<int>("-m"), rhfopt.thresh = hf.get<double>("-t");

        // gradient options
        gradhfopt.numerical = hf.get<std::vector<double>>("-g").at(0);
        gradhfopt.step = hf.get<std::vector<double>>("-g").at(1);

        // frequency options
        hesshfopt.numerical = hf.get<std::vector<double>>("-f").at(0);
        hesshfopt.step = hf.get<std::vector<double>>("-f").at(1);

        // optimization options
        opthfopt.thresh = hf.get<double>("-o");

        if (hf.is_subcommand_used("mp2")) {
            // printing options
            mp2print = mp2.get<std::vector<std::string>>("-p");

            // gradient options
            gradmpopt.numerical = mp2.get<std::vector<double>>("-g").at(0);
            gradmpopt.step = mp2.get<std::vector<double>>("-g").at(1);

            // frequency options
            gradmpopt.numerical = mp2.get<std::vector<double>>("-f").at(0);
            gradmpopt.step = mp2.get<std::vector<double>>("-g").at(1);

            // optimization options
            optmpopt.thresh = mp2.get<double>("-o");

        } else if (hf.is_subcommand_used("ci")) {
            // printing options
            ciprint = ci.get<std::vector<std::string>>("-p");
        }
    }

    // transform the print vectors to lowercase
    for (auto& element : mp2print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : hfprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : ciprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // calculate the integrals
    integrals();

    // distribute the calculations
    if (program.is_subcommand_used("hf")) hfrun();
}

void Distributor::hfrun() {
    // optimize the gradient
    if (hf.is_used("-o")) hfo();

    // print the RHF method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK METHOD\n" << std::string(104, '-') + "\n\n";
    std::printf("-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", rhfopt.maxiter, rhfopt.thresh, rhfopt.diis.start, rhfopt.diis.keep);

    // perform the Hartree-Fock calculation
    rhfres = HF(rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print the resulting matrices and energies
    if (CONTAINS(hfprint, "eps") || CONTAINS(print, "all")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(rhfres.eps) << std::endl;
    if (CONTAINS(hfprint, "c") || CONTAINS(print, "all")) std::cout << "\nCOEFFICIENT MATRIX\n" << rhfres.C << std::endl;
    if (CONTAINS(hfprint, "d") || CONTAINS(print, "all")) std::cout << "\nDENSITY MATRIX\n" << rhfres.D << std::endl;
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << rhfres.E << std::endl;

    // gradient and hessian frequency
    if (hf.is_used("-g")) hfg();
    if (hf.is_used("-f")) hff();

    // calculate the MP2 correlation;
    if (hf.is_subcommand_used("mp2")) mp2run();

    // calculate CI correlation
    if (hf.is_subcommand_used("ci")) {
        // print the CI method header and define J in MO basis
        std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;

        // transform the coulomb tensor
        std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
        if (CONTAINS(ciprint, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << Jmo;} std::cout << "\n";

        // do the calculation
        if (ci.get("-e") == "s") rcires = CI({rhfres}).cis(system, Jmo);
        if (ci.get("-e") == "d") rcires = CI({rhfres}).cid(system, Jmo);

        // print the result matrices
        if (CONTAINS(ciprint, "cih") || CONTAINS(print, "all")) std::cout << "\nCI HAMILTONIAN\n" << rcires.H << "\n";
        if (CONTAINS(ciprint, "cie") || CONTAINS(print, "all")) std::cout << "\nCI ENERGIES\n" << Matrix(rcires.eig) << "\n";
        if (CONTAINS(ciprint, "cic") || CONTAINS(print, "all")) std::cout << "\nCI EXPANSION COEFFICIENTS\n" << rcires.C << "\n";

        // print the gradient and norm
        std::cout << "\nCI CORRELATION ENERGY: " << rcires.Ecorr << std::endl;
        std::cout << "FINAL CI ENERGY: " << rhfres.E + rcires.Ecorr << std::endl;
    }
}

void Distributor::hff() {
    // print the hessian header
    if (hesshfopt.numerical) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n";

    // perform the frequency calculation
    auto[H, freq] = Hessian<HF>({rhfopt, rhfres, hf.get<std::vector<double>>("-f").at(0) == 1, hf.get<std::vector<double>>("-f").at(1)}).frequency(system);

    // print the hessian results
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());

    // print the frequency reslts
    std::cout << "\n" + std::string(104, '-') + "\nHARTREE-FOCK FREQUENCY ANALYSIS\n" << std::string(104, '-');
    std::cout << "\n\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::hfg() {
    // print the header
    if (gradhfopt.numerical) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK\n";
    std::cout << std::string(104, '-') + "\n\n"; 

    // calculate the gradient
    Matrix G = Gradient<HF>({rhfopt, rhfres, hf.get<std::vector<double>>("-g").at(0) == 1, hf.get<std::vector<double>>("-g").at(1)}).get(system);

    // print the results
    std::cout << G << "\n\nGRADIENT NORM: ";
    std::printf("%.2e\n", G.norm());
}

void Distributor::hfo() {
    // print the header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n\n";

    // perform the optimization
    Gradient<HF>::OptionsRestricted gradrhfopt = {rhfopt, rhfres, hf.get<std::vector<double>>("-g").at(0) == 1, hf.get<std::vector<double>>("-g").at(1)};
    system = Optimizer<HF>({rhfopt, gradrhfopt, hf.get<double>("-o")}).optimize(system);

    // print the results
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::mp2run() {
    // optimize the molecule with MP2 method
    if (mp2.is_used("-o")) mp2o();

    // print the MP2 correlation method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;

    // transform the coulomb tensor
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (CONTAINS(mp2print, "jmo") || CONTAINS(print, "all")) {std::cout << "\n" << Jmo;} std::cout << "\n";

    // do the MP2 calculation
    rmpres = MP({rhfres}).mp2(system, Jmo);

    // print the gradient and norm
    std::cout << "\nMP2 CORRELATION ENERGY: " << rmpres.Ecorr << std::endl << "FINAL ";
    std::cout << "MP2 ENERGY: " << rhfres.E + rmpres.Ecorr << std::endl;

    // calculate the MP2 nuclear gradient
    if (mp2.is_used("-g")) mp2g();
    if (mp2.is_used("-f")) mp2f();
}

void Distributor::mp2f() {
    // print the frequency calculation header
    if (hessmpopt.numerical) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL HESSIAN FOR RESTRICTED MP2 METHOD\n";
    std::cout << std::string(104, '-') + "\n\n";

    // perform the hessian calculation
    auto[H, freq] = Hessian<MP>({rhfopt, rhfres, mp2.get<std::vector<double>>("-f").at(0) == 1, mp2.get<std::vector<double>>("-f").at(1)}).frequency(system);

    // print the hessian resuls
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());

    // print the frequency analysis results
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 FREQUENCY ANALYSIS\n" << std::string(104, '-') + "\n";
    std::cout << "\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::mp2g() {
    // print the MP2 gradient method header and perform the calculation
    if (gradmpopt.numerical) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";

    // perform the calculation
    Matrix G = Gradient<MP>({rhfopt, rhfres, mp2.get<std::vector<double>>("-g").at(0) == 1, mp2.get<std::vector<double>>("-g").at(1)}).get(system);

    // print the gradient results
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::mp2o() {
    // print the MP2 optimization method header and optimize the molecule
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 OPTIMIZATION\n" << std::string(104, '-') << "\n\n";

    // perform the optimization
    Gradient<MP>::OptionsRestricted gradrmpopt = {rhfopt, rhfres, mp2.get<std::vector<double>>("-g").at(0) == 1, mp2.get<std::vector<double>>("-g").at(1)};
    system = Optimizer<MP>({rhfopt, gradrmpopt, mp2.get<double>("-o")}).optimize(system);

    // print the optimization results
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl;
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    if (CONTAINS(print, "s") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    if (CONTAINS(print, "t") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    if (CONTAINS(print, "v") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.V << std::endl;

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(print, "j") || CONTAINS(print, "all"))) {std::cout << "\n" << system.ints.J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if ((hf.is_used("-g") || hf.is_used("-o")) && !hf.get<std::vector<double>>("-g").at(0)) {
        // calculate the overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        if (CONTAINS(print, "ds") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        if (CONTAINS(print, "dt") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        if (CONTAINS(print, "dv") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        if (CONTAINS(print, "dj") || CONTAINS(print, "all")) {std::cout << "\n" << system.dints.dJ;} std::cout << "\n";
    }
}
