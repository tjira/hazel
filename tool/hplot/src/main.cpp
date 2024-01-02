#include "plot.h"

#include "argparse.hpp"

int main(int argc, char** argv) {
    // create the program
    argparse::ArgumentParser program("hazel", "0.1", argparse::default_arguments::none);

    // add arguments to program
    program.add_argument("-h", "--help").help("-- Print the help message.").default_value(false).implicit_value(true);
    program.add_argument("-mde").help("-- Energy from the MD simulation.").default_value(false).implicit_value(true);
    program.add_argument("-mdt").help("-- Temperature from the MD simulation.").default_value(false).implicit_value(true);
    program.add_argument("-hfc").help("-- Plot the convergence of the SCF loop.").default_value(false).implicit_value(true);
    program.add_argument("-pes").help("-- Plot the PES from a scan calculation.").default_value(false).implicit_value(true);

    // parse the arguments
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help message
    if (program.get<bool>("--help")) {std::cout << program; exit(EXIT_SUCCESS);}

    // read the piped input
    std::string input, line; while (getline(std::cin, line)) {input += line + "\n";} input.pop_back();

    // create the data container and print the input
    std::vector<std::vector<double>> data;
    std::cout << input << std::endl;

    if (program.get<bool>("-hfc")) {
        // create the line stringstream
        std::stringstream lss; lss << input;

        // loop over lines
        while (std::getline(lss, line)) {
            // if cursor in RHF block
            if (line == "RESTRICTED HARTREE-FOCK METHOD") {
                // skip lines
                for (int i = 0; i < 6; i++) std::getline(lss, line);

                // while iterations
                while (std::getline(lss, line) && !line.empty()) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // set precision
                    css << std::fixed; css.precision(14);

                    // exctract and append data
                    double iter, dE; css >> iter, css >> dE, css >> dE; data.push_back({iter, dE});
                }
            }
        }
    } else if (program.get<bool>("-mde")) {
        // create the line stringstream
        std::stringstream lss; lss << input;

        // loop over lines
        while (std::getline(lss, line)) {
            // if cursor in RHF block
            if (line == "MOLECULAR DYNAMICS") {
                // skip lines
                for (int i = 0; i < 6; i++) std::getline(lss, line);

                // while iterations
                while (std::getline(lss, line) && !line.empty()) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // set precision
                    css << std::fixed; css.precision(14);

                    // exctract and append data
                    double iter, time, E; css >> iter, css >> time, css >> E; data.push_back({time, E});
                }
            }
        }
    } else if (program.get<bool>("-mdt")) {
        // create the line stringstream
        std::stringstream lss; lss << input;

        // loop over lines
        while (std::getline(lss, line)) {
            // if cursor in RHF block
            if (line == "MOLECULAR DYNAMICS") {
                // skip lines
                for (int i = 0; i < 6; i++) std::getline(lss, line);

                // while iterations
                while (std::getline(lss, line) && !line.empty()) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // set precision
                    css << std::fixed; css.precision(14);

                    // exctract and append data
                    double iter, time, E, Ekin, T; css >> iter, css >> time, css >> E, css >> Ekin, css >> T; data.push_back({time, T});
                }
            }
        }
    } else if (program.get<bool>("-pes")) {
        // create the line stringstream
        std::stringstream lss; lss << input;

        // loop over lines
        while (std::getline(lss, line)) {
            // if cursor in RHF block
            if (line == "ENERGY SCAN") {
                // skip lines
                for (int i = 0; i < 3; i++) std::getline(lss, line);

                // while iterations
                while (std::getline(lss, line) && !line.empty()) {
                    // create the column stringstream
                    std::stringstream css; css << line;

                    // set precision
                    css << std::fixed; css.precision(14);

                    // exctract and append data
                    double iter, E; css >> iter, css >> E; data.push_back({iter, E});
                }
            }
        }
    }

    // create the gnuplot input
    Plot("plot").plot(data);
}
