#include "parser.h"

Parser::Parser(int argc, char** argv) : program("hazel", "0.1", argparse::default_arguments::none) {
    // reserve space for parsers

    // add main parsers
    parsers.reserve(6);
    parsers.insert({"ints", Parser("ints")}); program.add_subparser(parsers.at("ints").program);
    parsers.insert({"md", Parser("md")}); program.add_subparser(parsers.at("md").program);
    parsers.insert({"qd", Parser("qd")}); program.add_subparser(parsers.at("qd").program);
    parsers.insert({"rhf", Parser("rhf")}); program.add_subparser(parsers.at("rhf").program);
    parsers.insert({"scan", Parser("scan")}); program.add_subparser(parsers.at("scan").program);
    parsers.insert({"uhf", Parser("uhf")}); program.add_subparser(parsers.at("uhf").program);

    // add post RHF parsers
    at("rhf").parsers.reserve(6);
    at("rhf").parsers.insert({"cisd", Parser("cisd")}); at("rhf").program.add_subparser(at("rhf").at("cisd").program);
    at("rhf").parsers.insert({"mp2", Parser("mp2")}); at("rhf").program.add_subparser(at("rhf").at("mp2").program);
    at("rhf").parsers.insert({"cis", Parser("cis")}); at("rhf").program.add_subparser(at("rhf").at("cis").program);
    at("rhf").parsers.insert({"cid", Parser("cid")}); at("rhf").program.add_subparser(at("rhf").at("cid").program);
    at("rhf").parsers.insert({"fci", Parser("fci")}); at("rhf").program.add_subparser(at("rhf").at("fci").program);
    at("rhf").parsers.insert({"ci", Parser("ci")}); at("rhf").program.add_subparser(at("rhf").at("ci").program);

    // add parsers for MD
    at("md").parsers.reserve(2);
    at("md").parsers.insert({"rhf", Parser("rhf")}); at("md").program.add_subparser(at("md").at("rhf").program);
    at("md").parsers.insert({"uhf", Parser("uhf")}); at("md").program.add_subparser(at("md").at("uhf").program);

    // add MD post RHF parsers
    at("md").at("rhf").parsers.reserve(6);
    at("md").at("rhf").parsers.insert({"cisd", Parser("cisd")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cisd").program);
    at("md").at("rhf").parsers.insert({"mp2", Parser("mp2")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("mp2").program);
    at("md").at("rhf").parsers.insert({"cis", Parser("cis")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cis").program);
    at("md").at("rhf").parsers.insert({"cid", Parser("cid")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("cid").program);
    at("md").at("rhf").parsers.insert({"fci", Parser("fci")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("fci").program);
    at("md").at("rhf").parsers.insert({"ci", Parser("ci")}); at("md").at("rhf").program.add_subparser(at("md").at("rhf").at("ci").program);

    // add parsers for scan
    at("scan").parsers.reserve(2);
    at("scan").parsers.insert({"rhf", Parser("rhf")}); at("scan").program.add_subparser(at("scan").at("rhf").program);
    at("scan").parsers.insert({"uhf", Parser("uhf")}); at("scan").program.add_subparser(at("scan").at("uhf").program);

    // add scan post RHF parsers
    at("scan").at("rhf").parsers.reserve(6);
    at("scan").at("rhf").parsers.insert({"cisd", Parser("cisd")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cisd").program);
    at("scan").at("rhf").parsers.insert({"mp2", Parser("mp2")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("mp2").program);
    at("scan").at("rhf").parsers.insert({"cis", Parser("cis")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cis").program);
    at("scan").at("rhf").parsers.insert({"cid", Parser("cid")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("cid").program);
    at("scan").at("rhf").parsers.insert({"fci", Parser("fci")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("fci").program);
    at("scan").at("rhf").parsers.insert({"ci", Parser("ci")}); at("scan").at("rhf").program.add_subparser(at("scan").at("rhf").at("ci").program);

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
    program.at<argparse::ArgumentParser>("rhf").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-8).scan<'g', double>();

    // add arguments to the UHF argument parser
    program.at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-8).scan<'g', double>();

    // add arguments to the MP2 argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();

    // add arguments to the CI argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add arguments to the CIS argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add arguments to the CID argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add arguments to the CISD argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add arguments to the FCI argument parser
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add arguments to the MD argument parser
    program.at<argparse::ArgumentParser>("md").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").add_argument("-o", "--output").help("-- Output of the trajectory.").default_value("trajectory.xyz");

    // add arguments to the MD RHF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD UHF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("uhf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD MP2 argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD CI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD CIS argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD CID argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD CISD argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the MD FCI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add arguments to the scan argument parser
    program.at<argparse::ArgumentParser>("scan").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").add_argument("-o", "--output").help("-- Output of the PES energies.").default_value("pes.dat");

    // add arguments to the scan RHF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the scan UHF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("uhf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add arguments to the scan MP2 argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the scan CI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the scan CIS argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the scan CID argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the scan CISD argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the scan FCI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("rhf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add arguments to the QD argument parser
    program.at<argparse::ArgumentParser>("qd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("qd").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.1).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(1000).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-n", "--nstate").help("-- Number of states to consider.").default_value(3).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-o", "--output").help("-- Output of the wavefunction.").default_value("wavefunction.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-f", "--potfile").help("-- File with the PES.").default_value("pes.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-t", "--thresh").help("-- Threshold for conververgence in ITP loop.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("--no-real").help("-- Help message.").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& error) {
        std::cerr << error.what() << std::endl; exit(EXIT_FAILURE);
    }

    // print help of the base program if requested
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); exit(EXIT_SUCCESS);
    }

    // print help of all the subparsers
    // for (const auto& parser : parsers) {
    //     try {if (parser.get<bool>("-h")) {
    //         std::cout << parser.help().str(); exit(EXIT_SUCCESS);
    //     }} catch (const std::exception&) {}
    // }

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
