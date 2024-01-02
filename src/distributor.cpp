#include "distributor.h"

#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : parser(argc, argv), start(Timer::Now()) {}

Distributor::~Distributor() {
    std::cout << "\n"; Printer::Title(std::string("TOTAL EXECUTION TIME: ") + Timer::Format(Timer::Elapsed(start)).c_str());
}

void Distributor::run() {
    // initialize the system
    std::string basis = parser.get<std::string>("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
    std::ifstream stream(parser.get<std::string>("-f")); system = System(stream, basis, parser.get<int>("-c"), parser.get<int>("-s")); stream.close();

    // print the initial stuff
    Printer::Initial(parser, system);

    // print error if spin is not 1 for RHF
    if (parser.get<int>("-s") != 1) {
        if (parser.at("md").used("rhf")) {
            throw std::runtime_error("SPIN MUST BE 1 FOR RHF CALCULATIONS.");
        } else if (parser.at("scan").used("rhf")) {
            throw std::runtime_error("SPIN MUST BE 1 FOR RHF CALCULATIONS.");
        } else if (parser.used("rhf")) throw std::runtime_error("SPIN MUST BE 1 FOR RHF CALCULATIONS.");
    }

    // calculate the integrals if needed
    if (parser.used("ints") || parser.used("rhf") || parser.used("uhf")) integrals();

    // optimize the molecule
    if (parser.used("opt")) {
        if (parser.at("opt").used("rhf")) {
            if (parser.at("opt").at("rhf").used("mp2")) rmp2o();
            else if (parser.at("opt").at("rhf").used("ci")) rcio();
            else if (parser.at("opt").at("rhf").used("cis")) rcio();
            else if (parser.at("opt").at("rhf").used("cid")) rcio();
            else if (parser.at("opt").at("rhf").used("cisd")) rcio();
            else if (parser.at("opt").at("rhf").used("fci")) rcio();
            else rhfo();
        } else if (parser.at("opt").used("uhf")) uhfo();
        else throw std::runtime_error("NO METHOD FOR OPTIMIZATION SPECIFIED.");
    }

    // distribute the calculations
    if (parser.used("scan")) scan();
    if (parser.used("md")) dynamics();
    if (parser.used("qd")) qdyn();
    if (parser.used("rhf")) rhfrun();
    if (parser.used("uhf")) uhfrun();
}

void Distributor::rcirun(const HF::ResultsRestricted& rhfres) {
    // function to return the correct parser
    auto rciparser = [&]() -> Parser& {
        if (parser.at("rhf").used("ci")) return parser.at("rhf").at("ci");
        else if (parser.at("rhf").used("cid")) return parser.at("rhf").at("cid");
        else if (parser.at("rhf").used("cis")) return parser.at("rhf").at("cis");
        else if (parser.at("rhf").used("cisd")) return parser.at("rhf").at("cisd");
        else if (parser.at("rhf").used("fci")) return parser.at("rhf").at("fci");
        else throw std::runtime_error("INVALID CI METHOD");
    };

    // print the CI method header
    std::cout << "\n"; Printer::Title("RESTRICTED CONFIGURATION INTERACTION (" + Utility::ToUpper(rciparser().getName()) + ")");

    // define the matrices and results
    Matrix Hms; Tensor<4> Jms; CI::ResultsRestricted rcires;

    // transform the coulomb tensor
    std::cout << "\nCOULOMB INT IN MS BASIS: " << std::flush; TIME(Jms = Transform::CoulombSpin(system.ints.J, rhfres.C)) std::cout << " " << Eigen::MemTensor(Jms);
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "jms")) {std::cout << "\n" << Jms << "\n";} std::cout << "\n";
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-e"), "jms")) Eigen::Write("JMS.mat", Jms);
    
    // transform the hamiltonian
    std::cout << "HAMILTONIAN IN MS BASIS: " << std::flush; TIME(Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C)) std::cout << " " << Eigen::MemMatrix(Hms);
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "hms")) {std::cout << "\n" << Hms;} std::cout << "\n";
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-e"), "hms")) Eigen::Write("HMS.mat", Hms);

    // do the calculation
    if (parser.at("rhf").used("ci")) rcires = CI(rhfres).rci(system, parser.at("rhf").at("ci").get<std::vector<int>>("excitations"), Hms, Jms);
    if (parser.at("rhf").used("cisd")) rcires = CI(rhfres).rci(system, {1, 2}, Hms, Jms);
    if (parser.at("rhf").used("cis")) rcires = CI(rhfres).rci(system, {1}, Hms, Jms);
    if (parser.at("rhf").used("cid")) rcires = CI(rhfres).rci(system, {2}, Hms, Jms);
    if (parser.at("rhf").used("fci")) rcires = CI(rhfres).rci(system, {}, Hms, Jms);

    // print/save the result matrices
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "cih")) std::cout << "\nCI HAMILTONIAN\n" << rcires.H << "\n";
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "cie")) std::cout << "\nCI ENERGIES\n" << Matrix(rcires.eig) << "\n";
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "cic")) std::cout << "\nCI EXPANSION COEFFICIENTS\n" << rcires.C << "\n";
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-e"), "cih")) Eigen::Write("CIH.mat", rcires.H);
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-e"), "cie")) Eigen::Write("CIEPS.mat", Matrix(rcires.eig));
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-e"), "cic")) Eigen::Write("CIC.mat", rcires.C);

    // print the energy
    std::cout << "\nCI CORRELATION ENERGY: " << rcires.Ecorr << std::endl;
    std::cout << "FINAL CI ENERGY: " << rhfres.E + rcires.Ecorr << std::endl;

    // gradient and frequency
    if (rciparser().has("-g")) rcig(rhfres);
    if (rciparser().has("-f")) rcif(rhfres);
}

void Distributor::rcif(const HF::ResultsRestricted& rhfres) {
    // function to return the correct parser
    auto rciparser = [&]() -> Parser& {
        if (parser.at("rhf").used("ci")) return parser.at("rhf").at("ci");
        else if (parser.at("rhf").used("cid")) return parser.at("rhf").at("cid");
        else if (parser.at("rhf").used("cis")) return parser.at("rhf").at("cis");
        else if (parser.at("rhf").used("cisd")) return parser.at("rhf").at("cisd");
        else if (parser.at("rhf").used("fci")) return parser.at("rhf").at("fci");
        else throw std::runtime_error("INVALID CI METHOD");
    };

    // print the header for frequency calculation and define the gradient matrix
    if (rciparser().get<double>("-f")) {
        std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED CI (" + Utility::ToUpper(rciparser().getName()) + ")");
        std::printf("\n-- STEP: %.2e\n", rciparser().get<double>("-f"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED CI (" + Utility::ToUpper(rciparser().getName()) + ")");} Matrix H;

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb")); std::vector <int> excits;

    // choose the correct method
    if (parser.at("rhf").used("ci")) excits = rciparser().get<std::vector<int>>("excitations");
    if (parser.at("rhf").used("cisd")) excits = {1, 2};
    if (parser.at("rhf").used("cis")) excits = {1};
    if (parser.at("rhf").used("cid")) excits = {2};

    // perform the Hessian calculation
    if (rciparser().get<double>("-f")) H = Hessian(rciparser().get<double>("-f")).get(system, Lambda::ECI(rhfopt, excits, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RCI NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian(rciparser().get<double>("-f")).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED CI FREQUENCY ANALYSIS (" + Utility::ToUpper(rciparser().getName()) + ")");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rcig(const HF::ResultsRestricted& rhfres) {
    // function to return the correct parser
    auto rciparser = [&]() -> Parser& {
        if (parser.at("rhf").used("ci")) return parser.at("rhf").at("ci");
        else if (parser.at("rhf").used("cid")) return parser.at("rhf").at("cid");
        else if (parser.at("rhf").used("cis")) return parser.at("rhf").at("cis");
        else if (parser.at("rhf").used("cisd")) return parser.at("rhf").at("cisd");
        else if (parser.at("rhf").used("fci")) return parser.at("rhf").at("fci");
        else throw std::runtime_error("INVALID CI METHOD");
    };

    // print the header for gradient calculation and define the gradient matrix
    if (rciparser().get<double>("-g")) {
        std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED CI (" + Utility::ToUpper(rciparser().getName()) + ")");
        std::printf("\n-- STEP: %.2e\n", rciparser().get<double>("-g"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED CI (" + Utility::ToUpper(rciparser().getName()) + ")");} Matrix G; 

    // extract the RHF and RCI options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb")); std::vector <int> excits;

    // choose the correct method
    if (parser.at("rhf").used("ci")) excits = rciparser().get<std::vector<int>>("excitations");
    if (parser.at("rhf").used("cisd")) excits = {1, 2};
    if (parser.at("rhf").used("cis")) excits = {1};
    if (parser.at("rhf").used("cid")) excits = {2};

    // perform the CI gradient calculation
    if (rciparser().get<double>("-g")) G = Gradient(rciparser().get<double>("-g")).get(system, Lambda::ECI(rhfopt, excits, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RCI NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rcio() {
    // extract the RHF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("opt").at("rhf"), parser.get<bool>("--no-coulomb")); std::vector <int> excits;

    // function to return the correct parser
    auto rciparser = [&]() -> Parser& {
        if (parser.at("opt").at("rhf").used("ci")) return parser.at("opt").at("rhf").at("ci");
        else if (parser.at("opt").at("rhf").used("cid")) return parser.at("opt").at("rhf").at("cid");
        else if (parser.at("opt").at("rhf").used("cis")) return parser.at("opt").at("rhf").at("cis");
        else if (parser.at("opt").at("rhf").used("cisd")) return parser.at("opt").at("rhf").at("cisd");
        else if (parser.at("opt").at("rhf").used("fci")) return parser.at("opt").at("rhf").at("fci");
        else throw std::runtime_error("INVALID CI METHOD");
    };

    // print the RHF optimization header
    std::cout << std::endl; Printer::Title("RESTRICTED CI OPTIMIZATION (" + Utility::ToUpper(rciparser().getName()) + ")");

    // choose the correct method
    if (parser.at("opt").at("rhf").used("ci")) excits = rciparser().get<std::vector<int>>("excitations");
    if (parser.at("opt").at("rhf").used("cisd")) excits = {1, 2};
    if (parser.at("opt").at("rhf").used("cis")) excits = {1};
    if (parser.at("opt").at("rhf").used("cid")) excits = {2};

    // perform the optimization
    system = Optimizer(parser.at("opt").get<double>("-t")).optimize(system, Lambda::EGCI(rhfopt, excits, rciparser().get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(rciparser().get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::rhfrun() {
    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK METHOD");
    std::printf("\n-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", rhfopt.maxiter, rhfopt.thresh, rhfopt.diis.start, rhfopt.diis.keep);

    // perform the RHF calculation
    auto rhfres = HF(rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print/export the resulting RHF matrices and energies
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "eps")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(rhfres.eps) << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "c")) std::cout << "\nCOEFFICIENT MATRIX\n" << rhfres.C << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "d")) std::cout << "\nDENSITY MATRIX\n" << rhfres.D << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "eps")) Eigen::Write("EPS.mat", Matrix(rhfres.eps));
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "c")) Eigen::Write("C.mat", rhfres.C);
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "d")) Eigen::Write("D.mat", rhfres.D);

    // print population analysis
    if (Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "pop")) {
        std::cout << "\nMULLIKEN CHARGES\n" << Matrix(Wavetool::Mulliken(system, rhfres.D)) << std::endl;
    }

    // print the resulting energy
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << rhfres.E << std::endl;

    // gradient and frequency
    if (parser.at("rhf").has("-g")) rhfg(rhfres);
    if (parser.at("rhf").has("-f")) rhff(rhfres);

    // post RHF methods
    if (parser.at("rhf").used("cisd")) rcirun(rhfres);
    if (parser.at("rhf").used("mp2")) rmp2run(rhfres);
    if (parser.at("rhf").used("cis")) rcirun(rhfres);
    if (parser.at("rhf").used("cid")) rcirun(rhfres);
    if (parser.at("rhf").used("fci")) rcirun(rhfres);
    if (parser.at("rhf").used("ci")) rcirun(rhfres);
}

void Distributor::rhff(const HF::ResultsRestricted& rhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.at("rhf").get<double>("-f")) {
        std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED HARTREE-FOCK");
        std::printf("\n-- STEP: %.2e\n", parser.at("rhf").get<double>("-f"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED HARTREE-FOCK");} Matrix H; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.at("rhf").get<double>("-f")) H = Hessian(parser.at("rhf").get<double>("-f")).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RHF NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian(parser.at("rhf").get<double>("-f")).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rhfg(const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.at("rhf").get<double>("-g")) {
        std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK");
        std::printf("\n-- STEP: %.2e\n", parser.at("rhf").get<double>("-g"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK");} Matrix G; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.at("rhf").get<double>("-g")) G = Gradient(parser.at("rhf").get<double>("-g")).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else G = Gradient().get(system, rhfres) + Integral::dRepulsion(system);

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rhfo() {
    // print the RHF optimization header
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK OPTIMIZATION");

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("opt").at("rhf"), parser.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer(parser.at("opt").get<double>("-t")).optimize(system, Lambda::EGHF(rhfopt, parser.at("opt").at("rhf").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.at("opt").at("rhf").get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::uhfrun() {
    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("uhf"), parser.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK METHOD");
    std::printf("\n-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", uhfopt.maxiter, uhfopt.thresh, uhfopt.diis.start, uhfopt.diis.keep);

    // perform the RHF calculation
    auto uhfres = HF(uhfopt).uscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print/save the resulting RHF matrices and energies
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "epsa")) std::cout << "\nORBITAL ENERGIES (ALPHA)\n" << Matrix(uhfres.epsa) << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "epsb")) std::cout << "\nORBITAL ENERGIES (BETA)\n" << Matrix(uhfres.epsb) << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "ca")) std::cout << "\nCOEFFICIENT MATRIX (ALPHA)\n" << uhfres.Ca << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "cb")) std::cout << "\nCOEFFICIENT MATRIX (BETA)\n" << uhfres.Cb << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "da")) std::cout << "\nDENSITY MATRIX (ALPHA)\n" << uhfres.Da << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-p"), "db")) std::cout << "\nDENSITY MATRIX (BETA)\n" << uhfres.Db << std::endl;
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "epsa")) Eigen::Write("EPSA.mat", Matrix(uhfres.epsa));
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "epsb")) Eigen::Write("EPSB.mat", Matrix(uhfres.epsb));
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "ca")) Eigen::Write("CA.mat", uhfres.Ca);
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "cb")) Eigen::Write("CB.mat", uhfres.Cb);
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "da")) Eigen::Write("DA.mat", uhfres.Da);
    if (Utility::VectorContains<std::string>(parser.at("uhf").get<std::vector<std::string>>("-e"), "db")) Eigen::Write("DB.mat", uhfres.Db);
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << uhfres.E << std::endl;

    // gradient and frequency
    if (parser.at("uhf").has("-g")) uhfg(uhfres);
    if (parser.at("uhf").has("-f")) uhff(uhfres);
}

void Distributor::uhff(const HF::ResultsUnrestricted& uhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.at("uhf").get<double>("-f")) {
        std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR UNRESTRICTED HARTREE-FOCK");
        std::printf("\n-- STEP: %.2e\n", parser.at("uhf").get<double>("-f"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR UNRESTRICTED HARTREE-FOCK");} Matrix H; 

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("uhf"), parser.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.at("uhf").get<double>("-f")) H = Hessian(parser.at("uhf").get<double>("-f")).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR UHF NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian(parser.at("uhf").get<double>("-f")).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::uhfg(const HF::ResultsUnrestricted& uhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.at("uhf").get<double>("-g")) {
        std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK");
        std::printf("\n-- STEP: %.2e\n", parser.at("uhf").get<double>("-g"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK");} Matrix G; 

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("uhf"), parser.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.at("uhf").get<double>("-g")) G = Gradient(parser.at("uhf").get<double>("-g")).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::uhfo() {
    // print the UHF optimization header
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK OPTIMIZATION");

    // extract the HF options
    auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("opt").at("uhf"), parser.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer(parser.at("opt").get<double>("-t")).optimize(system, Lambda::EGHF(uhfopt, parser.at("opt").at("uhf").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.at("opt").at("uhf").get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::rmp2run(const HF::ResultsRestricted& rhfres) {
    // print the MP2 method header and declare JMO tensor
    std::cout << std::endl; Printer::Title("RESTRICTED MØLLER-PLESSET PERTRUBATION THEORY"); Tensor<4> Jmo;

    // transform the coulomb tensor to MO basis
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (Utility::VectorContains<std::string>(parser.at("rhf").at("mp2").get<std::vector<std::string>>("-p"), "jmo")) {std::cout << "\n" << Jmo;} std::cout << "\n";
    if (Utility::VectorContains<std::string>(parser.at("rhf").at("mp2").get<std::vector<std::string>>("-e"), "jmo")) Eigen::Write("JMO.mat", Jmo);

    // perform the MP2 calculation
    double Ecorr = MP(rhfres).rmp2(system, Jmo);

    // print the gradient and it's norm
    std::cout << "\nMP2 CORRELATION ENERGY: " << Ecorr << std::endl << "FINAL ";
    std::cout << "MP2 ENERGY: " << rhfres.E + Ecorr << std::endl;

    // gradient and frequency
    if (parser.at("rhf").at("mp2").has("-g")) rmp2g(rhfres);
    if (parser.at("rhf").at("mp2").has("-f")) rmp2f(rhfres);
}

void Distributor::rmp2f(const HF::ResultsRestricted& rhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.at("rhf").at("mp2").get<double>("-f")) {
        std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED MP2");
        std::printf("\n-- STEP: %.2e\n", parser.at("rhf").at("mp2").get<double>("-f"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED MP2");} Matrix H; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.at("rhf").at("mp2").get<double>("-f")) H = Hessian(parser.at("rhf").at("mp2").get<double>("-f")).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RMP2 NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian(parser.at("rhf").at("mp2").get<double>("-f")).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED MP2 FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rmp2g(const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.at("rhf").at("mp2").get<double>("-g")) {
        std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED MP2");
        std::printf("\n-- STEP: %.2e\n", parser.at("rhf").at("mp2").get<double>("-g"));
    } else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED MP2");} Matrix G; 

    // extract the RHF options
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("rhf"), parser.get<bool>("--no-coulomb"));

    // perform the MP2 gradient calculation
    if (parser.at("rhf").at("mp2").get<double>("-g")) G = Gradient(parser.at("rhf").at("mp2").get<double>("-g")).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RMP2 NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rmp2o() {
    // extract the HF options and print the MP2 optimization method header
    auto rhfopt = HF::OptionsRestricted::Load(parser.at("opt").at("rhf"), parser.get<bool>("--no-coulomb"));
    std::cout << std::endl; Printer::Title("RESTRICTED MP2 OPTIMIZATION");

    // perform the optimization
    system = Optimizer(parser.at("opt").get<double>("-t")).optimize(system, Lambda::EGMP2(rhfopt, parser.at("opt").at("rhf").at("mp2").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.at("opt").at("rhf").at("mp2").get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::scan() {
    // print the energy scan header
    std::cout << std::endl; Printer::Title("ENERGY SCAN"); std::printf("\nITER       Eel [Eh]           TIME    \n");

    // define the anonymous function for gradient
    std::function<double(System)> efunc; std::vector<double> energies;

    // get the energy and gradient function
    if (parser.at("scan").used("rhf")) {
        auto rhfopt = HF::OptionsRestricted::Load(parser.at("scan").at("rhf"), parser.get<bool>("--no-coulomb"));
        efunc = Lambda::EHF(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at("scan").at("rhf").used("mp2")) {
            efunc = Lambda::EMP2(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("scan").at("rhf").used("ci")) {
            efunc = Lambda::ECI(rhfopt, parser.at("scan").at("rhf").at("ci").get<std::vector<int>>("excitations"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("scan").at("rhf").used("cis")) {
            efunc = Lambda::ECI(rhfopt, {1}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("scan").at("rhf").used("cid")) {
            efunc = Lambda::ECI(rhfopt, {2}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("scan").at("rhf").used("cisd")) {
            efunc = Lambda::ECI(rhfopt, {1, 2}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("scan").at("rhf").used("fci")) {
            efunc = Lambda::ECI(rhfopt, {}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.at("scan").used("uhf")) {
        auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("scan").at("uhf"), parser.get<bool>("--no-coulomb"));
        efunc = Lambda::EHF(uhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
    }

    // count the number of geometries
    std::ifstream stream(parser.get<std::string>("-f")); std::string line; int lines;
    for (lines = 0; std::getline(stream, line); lines++);
    int geoms = lines / (system.atoms.size() + 2);
    stream.clear(), stream.seekg(0);

    // perform the scan
    for (int i = 0; i < geoms; i++) {
        // create the system
        std::string basis = parser.get<std::string>("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
        system = System(stream, basis, parser.get<int>("-c"), parser.get<int>("-s")); Timer::Timepoint start = Timer::Now();

        // start the timer and calculate the energy
        double E = efunc(system); energies.push_back(E);

        // print the energy
        std::printf("%4d %20.14f %s\n", i + 1, E, Timer::Format(Timer::Elapsed(start)).c_str());
    }
    
    // save the file
    std::ofstream file(parser.at("scan").get<std::string>("-o"));
    file << std::fixed << std::setprecision(14) << "# i              E\n";
    for (size_t i = 0; i < energies.size(); i++) {
        file << std::setw(5) << i << " " << std::setw(20) << energies.at(i) << "\n";
    }
}

void Distributor::dynamics() {
    // print the MD title
    std::cout << std::endl; Printer::Title("MOLECULAR DYNAMICS");

    // print MD header
    if (parser.at("md").has("--berendsen")) std::printf("\n-- THERMOSTAT: BERENDSEN, TEMP: %.2f, TAU: %.2f\n", parser.at("md").get<std::vector<double>>("--berendsen").at(0), parser.at("md").get<std::vector<double>>("--berendsen").at(1));
    else std::cout << std::endl;
    std::printf("-- ITERS: %d, TIMESTEP: %.2f\n", parser.at("md").get<int>("-i"), parser.at("md").get<double>("-s"));

    // define the anonymous function for gradient
    std::function<std::tuple<double, Matrix>(System&)> egfunc;

    // get the energy and gradient function
    if (parser.at("md").used("rhf")) {
        auto rhfopt = HF::OptionsRestricted::Load(parser.at("md").at("rhf"), parser.get<bool>("--no-coulomb"));
        egfunc = Lambda::EGHF(rhfopt, parser.at("md").at("rhf").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at("md").at("rhf").used("mp2")) {
            egfunc = Lambda::EGMP2(rhfopt, parser.at("md").at("rhf").at("mp2").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("md").at("rhf").used("ci")) {
            egfunc = Lambda::EGCI(rhfopt, parser.at("md").at("rhf").at("ci").get<std::vector<int>>("excitations"), parser.at("md").at("rhf").at("ci").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("md").at("rhf").used("cis")) {
            egfunc = Lambda::EGCI(rhfopt, {1}, parser.at("md").at("rhf").at("cis").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("md").at("rhf").used("cid")) {
            egfunc = Lambda::EGCI(rhfopt, {2}, parser.at("md").at("rhf").at("cid").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("md").at("rhf").used("cisd")) {
            egfunc = Lambda::EGCI(rhfopt, {1, 2}, parser.at("md").at("rhf").at("cisd").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at("md").at("rhf").used("fci")) {
            egfunc = Lambda::EGCI(rhfopt, {}, parser.at("md").at("rhf").at("fci").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.at("md").used("uhf")) {
        auto uhfopt = HF::OptionsUnrestricted::Load(parser.at("md").at("uhf"), parser.get<bool>("--no-coulomb"));
        egfunc = Lambda::EGHF(uhfopt, parser.at("md").at("uhf").get<double>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
    } else {
        throw std::runtime_error("NO METHOD SPECIFIED FOR MOLECULAR DYNAMICS");
    }

    // perform the dynamics
    MD(MD::Options::Load(parser.at("md"))).run(system, egfunc);
}

void Distributor::qdyn() {
    // print the dynamics title
    std::cout << std::endl; Printer::Title("QUANTUM DYNAMICS");

    // print QD header
    std::printf("\n-- ITERS: %d, TIMESTEP: %.2f, STATES: %d\n", parser.at("qd").get<int>("-i"), parser.at("qd").get<double>("-s"), parser.at("qd").get<int>("-n"));
    std::printf("-- THRESHOLD: %.2e\n\n", parser.at("qd").get<double>("-t"));

    // perform the dynamics
    QD::Results qdres = QD(QD::Options::Load(parser.at("qd"))).run(system);

    // print the energies
    if (parser.at("qd").get<bool>("--no-real")) std::cout << "\nIMAGINARY TIME PROPAGATION ENERGIES\n" << Matrix(qdres.energy) << std::endl;

    // save the wavefunction
    Utility::SaveWavefunction(parser.at("qd").get<std::string>("-o"), qdres.r, qdres.states, qdres.energy);
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n"; Printer::Title("INTEGRAL CALCULATION");

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.S);
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-e"), "s")) Eigen::Write("S.mat", system.ints.S);
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "s")) Eigen::Write("S.mat", system.ints.S);

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.T);
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-e"), "t")) Eigen::Write("T.mat", system.ints.T);
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "t")) Eigen::Write("T.mat", system.ints.T);

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.V);
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;
    if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-e"), "v")) Eigen::Write("V.mat", system.ints.V);
    if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "v")) Eigen::Write("V.mat", system.ints.V);

    // calculate the electron-electron repulsion integral
    if (!parser.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    if (!parser.get<bool>("--no-coulomb")) std::cout << " " << Eigen::MemTensor(system.ints.J);
    if (!parser.get<bool>("--no-coulomb") && parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;
    if (!parser.get<bool>("--no-coulomb") && parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;
    if (!parser.get<bool>("--no-coulomb") && parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-e"), "j")) Eigen::Write("J.mat", system.ints.J);
    if (!parser.get<bool>("--no-coulomb") && parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-e"), "j")) Eigen::Write("J.mat", system.ints.J);

    // print new line
    std::cout << "\n";

    // if derivatives of the integrals are needed
    if (parser.used("ints") || (parser.at("rhf").has("-g") && !parser.at("rhf").get<double>("-g"))) {
        // calculate the derivative of overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dS);
        if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;
        if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the derivative of kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dT);
        if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;
        if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the derivative of nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dV);
        if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;
        if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the derivative of electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dJ);
        if (parser.used("ints") && Utility::VectorContains<std::string>(parser.at("ints").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;
        if (parser.used("rhf") && Utility::VectorContains<std::string>(parser.at("rhf").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;

        // print new line
        std::cout << "\n";
    }
}
