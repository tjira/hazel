#include "distributor.h"

#define LOCALTIME [](){auto t = std::time(nullptr); auto tm = *std::localtime(&t); std::stringstream ss; ss << std::put_time(&tm, "%a %b %e %T %Y"); return ss.str();}()
#define CONTAINS(V, E) ([](std::vector<std::string> v, std::string e){return std::find(v.begin(), v.end(), e) != v.end();}(V, E))
#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}
#define ZERO Matrix::Zero(system.shells.nbf(), system.shells.nbf())
#define UNAME [](){utsname unr; uname(&unr); return unr;}()

Distributor::Distributor(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // initialize subcommands
    hf = argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none); mp2 = argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none);
    ci = argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none);

    // add positional arguments to the main argument parser
    program.add_argument("-b", "--basis").help("-- Basis set used to approximate atomic orbitals.").default_value(std::string("STO-3G"));
    program.add_argument("-c", "--charge").help("-- Charge of the system.").default_value(0).scan<'i', int>();
    program.add_argument("-f", "--file").help("-- Path to the system to use in the .xyz format.").default_value("molecule.xyz");
    program.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.add_argument("-n", "--nthread").help("-- Number of threads to use.").default_value(1).scan<'i', int>();
    program.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.add_argument("--no-center").help("-- Disable the molecule centering.").default_value(false).implicit_value(true);
    program.add_argument("--no-coulomb").help("-- Disable calculation of the coulomb tensor.").default_value(false).implicit_value(true);

    // add positional arguments to the HF argument parser
    hf.add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    hf.add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    hf.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    hf.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    hf.add_argument("-m", "--maxiter").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    hf.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    hf.add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    hf.add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    mp2.add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    mp2.add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    mp2.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    mp2.add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-8).scan<'g', double>();
    mp2.add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments tinput001.outo the CI argument parser
    ci.add_argument("-e", "--excitations").help("-- Excitations to consider.").default_value("d");
    ci.add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    ci.add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

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
    system = System(program.get("-f"), program.get("-b"), program.get<int>("-c"), 1);
    print = program.get<std::vector<std::string>>("-p");

    // check if unrestricted calculation needed
    if (system.charge % 2) throw std::runtime_error("SPIN UNRESTRICTED CALCULATIONS ARE NOT SUPPORTED YET");

    // print the title with number of threads
    std::cout << "QUANTUM HAZEL" << std::endl;

    // print the general info block
    std::cout << "\n" + std::string(104, '-') + "\nGENERAL INFO\n" << std::string(104, '-') + "\n\n";
    std::printf("COMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, LOCALTIME.c_str());
    std::printf("\nCOMPILER VERSION AND FLAGS: MACHINE %s, GCC %d.%d.%d (%s)", UNAME.machine, __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
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
    if (program.is_subcommand_used("hf")) {
        hfprint = hf.get<std::vector<std::string>>("-p");
        if (hf.is_subcommand_used("mp2")) {
            mp2print = mp2.get<std::vector<std::string>>("-p");
        } else if (hf.is_subcommand_used("ci")) {
            ciprint = ci.get<std::vector<std::string>>("-p");
        }
    }

    // extract method options
    if (program.is_subcommand_used("hf")) {
        rhfopt.diis = {hf.get<std::vector<int>>("-d").at(0), hf.get<std::vector<int>>("-d").at(1)};
        rhfopt.maxiter = hf.get<int>("-m"), rhfopt.thresh = hf.get<double>("-t");
        rhfopt.nocoulomb = program.get<bool>("--no-coulomb");
    }

    // transform the print vectors to lowercase
    for (auto& element : mp2print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : hfprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : ciprint) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});
    for (auto& element : print) std::transform(element.begin(), element.end(), element.begin(), [](auto c){return std::tolower(c);});

    // calculate the integrals
    integrals();

    // optimize the molecule
    if (program.is_subcommand_used("hf")) {
        if (hf.is_used("-o")) rhfo();
        else if (hf.is_subcommand_used("mp2")) {
            if (mp2.is_used("-o")) rmp2o();
        }
    }

    // distribute the calculations
    if (program.is_subcommand_used("hf")) rhfrun();
}

void Distributor::rcirun() {
    // print the CI method header and define J in MO basis
    std::cout << "\n" + std::string(104, '-') + "\nCI CORRELATION ENERGY\n" << std::string(104, '-') + "\n"; Tensor<4> Jmo;

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
    rhfres = HF(rhfopt).rscf(system, ZERO);

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
    std::cout << std::string(104, '-') + "\n\n";

    // define the energy function for numerical Hessian
    auto efunc = [this](System system) {
        return HF(rhfopt).rscf(system.clearints(), rhfres.D, false).E;
    };

    // perform the Hessian calculation
    Matrix H = Hessian({hf.get<std::vector<double>>("-f").at(1)}).get(system, efunc);

    // print the Hessian and it's norm
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());

    // perform the frequency calculation based on the calculated Hessian
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

    // define the energy function for numerical gradient
    auto efunc = [this](System system) {
        return HF(rhfopt).rscf(system.clearints(), rhfres.D, false).E;
    };

    // calculate the numerical or analytical gradient
    if (hf.get<std::vector<double>>("-g").at(0)) G = Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, efunc);
    else G = Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, rhfres);

    // print the gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rhfo() {
    // print the RHF optimization header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED HARTREE-FOCK OPTIMIZATION\n" << std::string(104, '-') + "\n\n";

    // define the energy function formnumerical gradient
    auto efunc = [this](System system) {
        return HF(rhfopt).rscf(system.clearints(), Matrix::Zero(system.shells.nbf(), system.shells.nbf()), false).E;
    };

    // define energy-gradient function
    auto egfunc = [this, efunc](System& system) {
        // delete the calculated integrals and recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), ZERO, false);

        // calculate the numerical or analytical gradient
        if (hf.get<std::vector<double>>("-g").at(0)) return std::tuple{rhfres.E, Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, efunc, false)};
        else return std::tuple{rhfres.E, Gradient({hf.get<std::vector<double>>("-g").at(1)}).get(system, rhfres, false)};
    };

    // perform the optimization
    system = Optimizer({hf.get<double>("-o")}).optimize(system, egfunc);

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl; 
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
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
    std::cout << std::string(104, '-') + "\n\n";

    // define the energy function for numerical Hessian
    auto efunc = [this](System system) {
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), this->rhfres.D, false);
        return rhfres.E + MP({rhfres}).rmp2(system, Tensor<4>(), false);
    };

    // perform the Hessian calculation
    Matrix H = Hessian({mp2.get<std::vector<double>>("-f").at(1)}).get(system, efunc);

    // print the Hessian and it's norm
    std::cout << "NUCLEAR HESSIAN\n" << H << "\n\nHESSIAN NORM: "; std::printf("%.2e\n", H.norm());

    // perform the frequency calculation based on the calculated Hessian
    Vector freq = Hessian({mp2.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the MP2 frequencies
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 FREQUENCY ANALYSIS\n" << std::string(104, '-') + "\n";
    std::cout << "\nVIBRATIONAL FREQUENCIES\n" << Matrix(freq) << std::endl;
}

void Distributor::rmp2g() {
    // print the MP2 gradient method header
    if (hf.get<std::vector<double>>("-g").at(0)) std::cout << "\n" + std::string(104, '-') + "\nNUMERICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";
    else std::cout << "\n" + std::string(104, '-') + "\nANALYTICAL GRADIENT FOR RESTRICTED MP2 METHOD\n" << std::string(104, '-') << "\n\n";

    // define the energy function for numerical gradient
    auto efunc = [this](System system) {
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), this->rhfres.D, false);
        return rhfres.E + MP({rhfres}).rmp2(system, Tensor<4>(), false);
    };

    // perform the MP2 gradient calculation
    Matrix G = Gradient({mp2.get<std::vector<double>>("-g").at(1)}).get(system, efunc);

    // print the MP2 gradient and it's norm
    std::cout << G << "\n\nGRADIENT NORM: "; std::printf("%.2e\n", G.norm());
}

void Distributor::rmp2o() {
    // print the MP2 optimization method header
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED MP2 OPTIMIZATION\n" << std::string(104, '-') << "\n\n";

    // define the energy function for numerical gradient
    auto efunc = [this](System system) {
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), ZERO, false);
        return rhfres.E + MP({rhfres}).rmp2(system, Tensor<4>(), false);
    };

    // define the energy-gradient function
    auto egfunc = [this, efunc](System& system) {
        // delete the calculated integrals and recalculate Hartree-Fock
        HF::ResultsRestricted rhfres = HF(rhfopt).rscf(system.clearints(), ZERO, false);

        // calculate the MP2 gradient and energy
        Matrix G = Gradient({mp2.get<std::vector<double>>("-g").at(1)}).get(system, efunc, false);
        double Ecorr = MP({rhfres}).rmp2(system, Tensor<4>(), false);

        // return the tuple containing energy and gradient
        return std::tuple{rhfres.E + Ecorr, G};
    };

    // perform the optimization
    system = Optimizer({mp2.get<double>("-o")}).optimize(system, egfunc);

    // print the optimized coordinates and distances
    std::cout << "\nOPTIMIZED SYSTEM COORDINATES\n" << system.coords << std::endl;
    if (CONTAINS(print, "dist") || CONTAINS(print, "all")) std::cout << "\nOPTIMIZED DISTANCE MATRIX\n" << system.dists << std::endl;
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n" + std::string(104, '-') + "\nINTEGRAL CALCULATION\n";
    std::cout << std::string(104, '-') << std::endl;

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.S);
    if (CONTAINS(print, "s") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.S << std::endl;

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.T);
    if (CONTAINS(print, "t") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.T << std::endl;

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.V);
    if (CONTAINS(print, "v") || CONTAINS(print, "all")) std::cout << "\n" << system.ints.V << std::endl;

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    std::cout << " " << Eigen::MemTensor(system.ints.J);
    if (!program.get<bool>("--no-coulomb") && (CONTAINS(print, "j") || CONTAINS(print, "all"))) {std::cout << "\n" << system.ints.J;} std::cout << "\n";

    // if derivatives of the integrals are needed
    if ((hf.is_used("-g") || hf.is_used("-o")) && !hf.get<std::vector<double>>("-g").at(0)) {
        // calculate the derivative of overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dS);
        if (CONTAINS(print, "ds") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the derivative of kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dT);
        if (CONTAINS(print, "dt") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the derivative of nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dV);
        if (CONTAINS(print, "dv") || CONTAINS(print, "all")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the derivative of electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dJ);
        if (CONTAINS(print, "dj") || CONTAINS(print, "all")) {std::cout << "\n" << system.dints.dJ;} std::cout << "\n";
    }
}
