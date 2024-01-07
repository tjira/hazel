#include "parser.h"

Parser::Parser(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none) {
    // add main parsers
    parsers.reserve(9);
    parsers.insert({"bagel", std::make_shared<Parser>(Parser("bagel"))}); program.add_subparser(at("bagel").program);
    parsers.insert({"ints", std::make_shared<Parser>(Parser("ints"))}); program.add_subparser(at("ints").program);
    parsers.insert({"md", std::make_shared<Parser>(Parser("md"))}); program.add_subparser(at("md").program);
    parsers.insert({"opt", std::make_shared<Parser>(Parser("opt"))}); program.add_subparser(at("opt").program);
    parsers.insert({"orca", std::make_shared<Parser>(Parser("orca"))}); program.add_subparser(at("orca").program);
    parsers.insert({"qd", std::make_shared<Parser>(Parser("qd"))}); program.add_subparser(at("qd").program);
    parsers.insert({"rhf", std::make_shared<Parser>(Parser("rhf"))}); program.add_subparser(at("rhf").program);
    parsers.insert({"scan", std::make_shared<Parser>(Parser("scan"))}); program.add_subparser(at("scan").program);
    parsers.insert({"uhf", std::make_shared<Parser>(Parser("uhf"))}); program.add_subparser(at("uhf").program);

    // add post RHF parsers
    at("rhf").parsers.reserve(6);
    at("rhf").parsers.insert({"cisd", std::make_shared<Parser>(Parser("cisd"))}); at("rhf").program.add_subparser(at("rhf").at("cisd").program);
    at("rhf").parsers.insert({"mp2", std::make_shared<Parser>(Parser("mp2"))}); at("rhf").program.add_subparser(at("rhf").at("mp2").program);
    at("rhf").parsers.insert({"cis", std::make_shared<Parser>(Parser("cis"))}); at("rhf").program.add_subparser(at("rhf").at("cis").program);
    at("rhf").parsers.insert({"cid", std::make_shared<Parser>(Parser("cid"))}); at("rhf").program.add_subparser(at("rhf").at("cid").program);
    at("rhf").parsers.insert({"fci", std::make_shared<Parser>(Parser("fci"))}); at("rhf").program.add_subparser(at("rhf").at("fci").program);
    at("rhf").parsers.insert({"ci", std::make_shared<Parser>(Parser("ci"))}); at("rhf").program.add_subparser(at("rhf").at("ci").program);

    // add parsers for OPT
    at("opt").parsers.reserve(4);
    at("opt").parsers.insert({"bagel", std::make_shared<Parser>(Parser("bagel"))}); at("opt").program.add_subparser(at("opt").at("bagel").program);
    at("opt").parsers.insert({"orca", std::make_shared<Parser>(Parser("orca"))}); at("opt").program.add_subparser(at("opt").at("orca").program);
    at("opt").parsers.insert({"rhf", std::make_shared<Parser>(Parser("rhf"))}); at("opt").program.add_subparser(at("opt").at("rhf").program);
    at("opt").parsers.insert({"uhf", std::make_shared<Parser>(Parser("uhf"))}); at("opt").program.add_subparser(at("opt").at("uhf").program);

    // add OPT post RHF parsers
    at("opt").at("rhf").parsers.reserve(6);
    at("opt").at("rhf").parsers.insert({"cisd", std::make_shared<Parser>(Parser("cisd"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("cisd").program);
    at("opt").at("rhf").parsers.insert({"mp2", std::make_shared<Parser>(Parser("mp2"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("mp2").program);
    at("opt").at("rhf").parsers.insert({"cis", std::make_shared<Parser>(Parser("cis"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("cis").program);
    at("opt").at("rhf").parsers.insert({"cid", std::make_shared<Parser>(Parser("cid"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("cid").program);
    at("opt").at("rhf").parsers.insert({"fci", std::make_shared<Parser>(Parser("fci"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("fci").program);
    at("opt").at("rhf").parsers.insert({"ci", std::make_shared<Parser>(Parser("ci"))}); at("opt").at("rhf").program.add_subparser(at("opt").at("rhf").at("ci").program);

    // add parsers for MD
    at("md").parsers.reserve(4);
    at("md").parsers.insert({"bagel", std::make_shared<Parser>(Parser("bagel"))}); at("md").program.add_subparser(at("md").at("bagel").program);
    at("md").parsers.insert({"orca", std::make_shared<Parser>(Parser("orca"))}); at("md").program.add_subparser(at("md").at("orca").program);
    at("md").parsers.insert({"rhf", std::make_shared<Parser>(Parser("rhf"))}); at("md").program.add_subparser(at("md").at("rhf").program);
    at("md").parsers.insert({"uhf", std::make_shared<Parser>(Parser("uhf"))}); at("md").program.add_subparser(at("md").at("uhf").program);

    // add MD post RHF parsers
    at("md").at("rhf").parsers.reserve(6);
    at("md").at("rhf").parsers.insert({"cisd", std::make_shared<Parser>(Parser("cisd"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cisd").program);
    at("md").at("rhf").parsers.insert({"mp2", std::make_shared<Parser>(Parser("mp2"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("mp2").program);
    at("md").at("rhf").parsers.insert({"cis", std::make_shared<Parser>(Parser("cis"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cis").program);
    at("md").at("rhf").parsers.insert({"cid", std::make_shared<Parser>(Parser("cid"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cid").program);
    at("md").at("rhf").parsers.insert({"fci", std::make_shared<Parser>(Parser("fci"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("fci").program);
    at("md").at("rhf").parsers.insert({"ci", std::make_shared<Parser>(Parser("ci"))}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("ci").program);

    // add parsers for SCAN
    at("scan").parsers.reserve(4);
    at("scan").parsers.insert({"bagel", std::make_shared<Parser>(Parser("bagel"))}); at("scan").program.add_subparser(at("scan").at("bagel").program);
    at("scan").parsers.insert({"orca", std::make_shared<Parser>(Parser("orca"))}); at("scan").program.add_subparser(at("scan").at("orca").program);
    at("scan").parsers.insert({"rhf", std::make_shared<Parser>(Parser("rhf"))}); at("scan").program.add_subparser(at("scan").at("rhf").program);
    at("scan").parsers.insert({"uhf", std::make_shared<Parser>(Parser("uhf"))}); at("scan").program.add_subparser(at("scan").at("uhf").program);

    // add scan post RHF parsers
    at("scan").at("rhf").parsers.reserve(6);
    at("scan").at("rhf").parsers.insert({"cisd", std::make_shared<Parser>(Parser("cisd"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cisd").program);
    at("scan").at("rhf").parsers.insert({"mp2", std::make_shared<Parser>(Parser("mp2"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("mp2").program);
    at("scan").at("rhf").parsers.insert({"cis", std::make_shared<Parser>(Parser("cis"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cis").program);
    at("scan").at("rhf").parsers.insert({"cid", std::make_shared<Parser>(Parser("cid"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cid").program);
    at("scan").at("rhf").parsers.insert({"fci", std::make_shared<Parser>(Parser("fci"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("fci").program);
    at("scan").at("rhf").parsers.insert({"ci", std::make_shared<Parser>(Parser("ci"))}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("ci").program);

    // add arguments to the main argument parser
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

    // add arguments to the integral argument parser
    program.at<argparse::ArgumentParser>("ints").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("ints").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("ints").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the RHF argument parser
    program.at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the UHF argument parser
    program.at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the MP2 argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the CI argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the CIS argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the CID argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the CISD argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the FCI argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to OPT argument parser
    program.at<argparse::ArgumentParser>("opt").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").add_argument("-t", "--threshold").help("-- Gradient threshold.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT RHF argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT UHF argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("uhf").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT MP2 argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT CI argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT CIS argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT CID argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT CISD argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT FCI argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the OPT BAGEL argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("bagel").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("bagel").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("bagel").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the OPT ORCA argument parser
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("orca").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("orca").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("opt").at<argparse::ArgumentParser>("orca").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the MD argument parser
    program.at<argparse::ArgumentParser>("md").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").add_argument("-e", "--excitation").help("-- Initial state of the molecule.").default_value(1).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").add_argument("-s", "--step").help("-- Time step in atomic units.").default_value(1.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").add_argument("-o", "--output").help("-- Output of the trajectory.").default_value("trajectory.xyz");
    program.at<argparse::ArgumentParser>("md").add_argument("--berendsen").help("-- Enable Berendsen thermostat with temperature and tau parameter.").default_value(std::vector<double>{0.0, 1.0}).nargs(2).scan<'g', double>();

    // add arguments to the MD RHF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();

    // add arguments to the MD UHF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD MP2 argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD CI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD CIS argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD CID argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD CISD argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD FCI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(1e-5).scan<'g', double>();

    // add arguments to the MD BAGEL argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("bagel").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("bagel").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("bagel").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the MD ORCA argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("orca").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("orca").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("orca").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the SCAN argument parser
    program.at<argparse::ArgumentParser>("scan").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").add_argument("-o", "--output").help("-- Output of the PES energies.").default_value("pes.dat");
    program.at<argparse::ArgumentParser>("scan").add_argument("-n", "--nstate").help("-- Number of scanned states.").default_value(1).scan<'i', int>();

    // add arguments to the SCAN RHF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the SCAN UHF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the SCAN MP2 argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN CI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Exctitations to use for the CI calculation.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN CIS argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN CID argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN CISD argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN FCI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the SCAN BAGEL argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("bagel").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("bagel").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the SCAN ORCA argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("orca").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("orca").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the QD argument parser
    program.at<argparse::ArgumentParser>("qd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("qd").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.1).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(1000).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-n", "--nstate").help("-- Number of states to consider.").default_value(3).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-o", "--output").help("-- Output of the wavefunction.").default_value("wavefunction.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-f", "--potfile").help("-- File with the PES.").default_value("pes.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-t", "--thresh").help("-- Threshold for conververgence in ITP loop.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("--no-real").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the BAGEL argument parser
    program.at<argparse::ArgumentParser>("bagel").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("bagel").add_argument("-g", "--gradient").help("-- States to calculate the gradient for.").nargs(argparse::nargs_pattern::any).default_value(std::vector<int>{0}).scan<'i', int>();
    program.at<argparse::ArgumentParser>("bagel").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("bagel").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    // add arguments to the ORCA argument parser
    program.at<argparse::ArgumentParser>("orca").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("orca").add_argument("-g", "--gradient").help("-- Step size for gradient calculation or 0 for analytical gradient.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("orca").add_argument("-f", "--frequency").help("-- Step size for frequency calculation or 0 for analytical hessian.").default_value(0.0).scan<'g', double>();
    program.at<argparse::ArgumentParser>("orca").add_argument("-m", "--method").help("-- Method for ORCA calculation.").default_value("hf");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    } help();

    // print all tha available bases if requested
    if (program.get<bool>("--show-bases")) {
        for (const auto& file : std::filesystem::directory_iterator(std::string(DATADIR) + "/basis")) {
            if (file.path().filename() != "basis.sh") std::cout << file.path().filename().replace_extension() << std::endl;
        }
        exit(EXIT_SUCCESS);
    }

    // set the path to the basis functions
    if (!std::filesystem::is_directory(std::string(DATADIR) + "/basis")) {
        auto path = std::filesystem::weakly_canonical(std::filesystem::path(argv[0])).parent_path();
        #ifdef _WIN32
        _putenv_s("LIBINT_DATA_PATH", path.string().c_str());
        #else
        setenv("LIBINT_DATA_PATH", path.c_str(), true);
        #endif
    }
    
    // set number of threads to use and cout flags
    std::cout << std::fixed << std::setprecision(14);
    nthread = program.get<int>("--nthread");
}

void Parser::help() const {
    if (has("-h")) {std::cout << program.help().str(); exit(EXIT_SUCCESS);}
    for (const auto& parser : parsers) parser.second->help();
}
