#include "argparse.hpp"

int main(int argc, char** argv) {
    // create the program
    argparse::ArgumentParser program("hazel", "0.1", argparse::default_arguments::none);

    // add arguments to program
    program.add_argument("-h", "--help").help("-- Print the help message.").default_value(false).implicit_value(true);
    program.add_argument("--mdenergy").help("-- Energy from the MD simulation.").default_value(false).implicit_value(true);
    program.add_argument("--hfconv").help("-- Plot the convergence of the SCF loop.").default_value(false).implicit_value(true);

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

    // print the input
    std::cout << input << std::endl;
}
