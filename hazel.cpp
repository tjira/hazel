#include "include/hartreefock.h"
#include "include/molecule.h"
#include <boost/json.hpp>
#include <boost/program_options.hpp>
#include <filesystem>

namespace js = boost::json;
namespace po = boost::program_options;

int main(int argc, char** argv) {
    po::options_description desc("options");
    po::positional_options_description pos;
    po::variables_map vm;

    desc.add_options()
        ("input", po::value<std::string>(), "input file")
        ("version,v", "print version string")
        ("help,h", "produce help message")
    ;pos.add("input", 1);

    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm); po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl; return EXIT_SUCCESS;
    }

    std::ifstream file(vm["input"].as<std::string>());
    std::stringstream buffer; buffer << file.rdbuf();
    js::value input = js::parse(buffer.str());

    std::string path = std::filesystem::absolute(vm["input"].as<std::string>()).parent_path();
    std::string molfile = js::value_to<std::string>(input.at("molecule"));
    std::string basis = js::value_to<std::string>(input.at("basis"));
    
    HartreeFockOptions opt = {
        js::value_to<double>(input.at("damp")),
        js::value_to<double>(input.at("thresh")),
        js::value_to<int>(input.at("maxiter"))
    };
    Molecule molecule(path + "/" + molfile, basis);
    Printer::printTitle();
    HartreeFock hfock(opt);
    
    libint2::initialize();
    auto result = hfock.scf(molecule);
    libint2::finalize();
}
