#include "printer.h"

void Printer::Initial(const Parser& program, System& system) {
    // print the title
    std::cout << "QUANTUM HAZEL\n\n";

    // print the general info block
    Title("GENERAL INFO");
    std::printf("\nCOMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, Timer::Local().c_str());
    std::printf("\nCOMPILER: %s, GCC %d.%d.%d (%s)", OS, __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);

    // print the external libraries
    #ifndef _WIN32
    // get the ORCA and BAGEL versions
    bool bagel = !std::system("which BAGEL > /dev/null 2>&1");
    bool orca = !std::system("which orca > /dev/null 2>&1");

    // check if any of the two is available
    if (orca || bagel) {
        // print the external libraries block
        std::printf("\nEXTERNAL:");

        // print the ORCA version
        if (orca) {
            // define the command from which to get the ORCA version
            std::unique_ptr<FILE, decltype(&pclose)> pipe(popen("orca - 2>&1", "r"), pclose);

            // check for success
            if (!pipe) throw std::runtime_error("ORCA EXECUTION FAILED");

            // read the output 
            std::array<char, 128> buffer; std::string output;
            while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
                output += std::string(buffer.data());
            }

            // extracto orca version using regex
            std::smatch match; std::regex_search(output, match, std::regex("Program Version ([0-9]+)\\.([0-9]+)\\.([0-9]+)"));

            // print the version
            std::printf(" ORCA %s", match.str(0).substr(16, 5).c_str());
        }

        // print the BAGEL version
        if (bagel) {
            std::cout << (orca ? ", " : " ") << "BAGEL";
        }

        // print the newline
        std::cout << std::endl;
    } else {
        std::cout << std::endl;
    }
    #else
    std::cout << std::endl;
    #endif

    // print the OpenMP info
    #if defined(_OPENMP)
    std::cout << std::endl; Title(std::string("OPENMP ") + std::to_string(_OPENMP));
    std::printf("\nAVAILABLE CORES: %d\nUSED THREADS: %d\n", std::thread::hardware_concurrency(), program.get<int>("-n"));
    #endif

    if (!program.used("qd")) {
        // empty line
        std::cout << std::endl;

        // print the system block
        Title(std::string("SYSTEM SPECIFICATION (") + Utility::ToUpper(program.get<std::string>("-b")) + std::string(")"));
        std::printf("\n-- ATOMS: %d, ELECTRONS: %d, NBF: %d\n", (int)system.atoms.size(), system.electrons, (int)system.shells.nbf());
        std::printf("-- CHARGE: %d, MULTIPLICITY: %d\n", system.charge, system.multi);
        Printer::Mat("\nSYSTEM COORDINATES", system.coords);

        // center the molecule if not disabled
        if (!program.get<bool>("--no-center")) {
            Matrix dir(system.atoms.size(), 3); dir.rowwise() -= system.coords.colwise().sum() / system.atoms.size();
            system.move(dir * A2BOHR); Printer::Mat("\nCENTERED SYSTEM COORDINATES", system.coords);
        }

        // print the distances if requested
        if (Utility::VectorContains<std::string>(program.get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nDISTANCE MATRIX", system.dists);
    }
}

void Printer::Mat(const std::string& title, const Matrix& A) {
    std::cout << title << "\n" << A << std::endl; 
}

void Printer::Title(const std::string& title) {
    std::cout << std::string(124, '-') + "\n" << title << "\n" << std::string(124, '-') << std::endl;
}
