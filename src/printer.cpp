#include "printer.h"

void Printer::Initial(const Parser& program, System& system) {
    // print the title
    std::cout << "QUANTUM HAZEL\n\n";

    // print the general info block
    Title("GENERAL INFO");
    std::printf("\nCOMPILATION TIMESTAMP: %s\nEXECUTION TIMESTAMP: %s\n", __TIMESTAMP__, Timer::Local().c_str());
    std::printf("\nCOMPILER: %s, GCC %d.%d.%d (%s)", OS, __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__, CXXFLAGS);
    std::printf("\nLIBRARIES: EIGEN %d.%d.%d, LIBINT %d.%d.%d\n", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, LIBINT_MAJOR_VERSION, LIBINT_MINOR_VERSION, LIBINT_MICRO_VERSION);
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
    std::cout << std::string(104, '-') + "\n" << title << "\n" << std::string(104, '-') << std::endl;
}
