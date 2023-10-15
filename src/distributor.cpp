#include "distributor.h"

#define TIME(W) {Timer::Timepoint start = Timer::Now(); W; std::cout << Timer::Format(Timer::Elapsed(start)) << std::flush;}

Distributor::Distributor(int argc, char** argv) : parsers(25), program("hazel", "0.1", argparse::default_arguments::none), start(Timer::Now()) {
    // add level 1 parsers
    parsers.push_back(argparse::ArgumentParser("ints", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("scan", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("qd", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("md", "0.1", argparse::default_arguments::none)); program.add_subparser(parsers.at(parsers.size() - 1));

    // add HF parsers
    parsers.push_back(argparse::ArgumentParser("cisd", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cis", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cid", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("fci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));

    // add side parsers
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("hf", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").add_subparser(parsers.at(parsers.size() - 1));

    // add scan parsers
    parsers.push_back(argparse::ArgumentParser("cisd", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cis", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cid", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("fci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));

    // add MD parsers
    parsers.push_back(argparse::ArgumentParser("cisd", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("mp2", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cis", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("cid", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("fci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    parsers.push_back(argparse::ArgumentParser("ci", "0.1", argparse::default_arguments::none)); program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_subparser(parsers.at(parsers.size() - 1));
    
    // add positional arguments to the main argument parser
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

    // add positional arguments to the integral argument parser
    program.at<argparse::ArgumentParser>("ints").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("ints").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("ints").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the HF argument parser
    program.at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").add_argument("-o", "--optimize").help("-- Optimization with gradient threshold.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-8).scan<'g', double>();

    // add positional arguments to the MP2 argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-p", "--print").help("-- Output printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();

    // add positional arguments to the CI argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add positional arguments to the CIS argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add positional arguments to the CID argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add positional arguments to the CISD argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add positional arguments to the FCI argument parser
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-f", "--frequency").help("-- Analytical (0) or numerical (1) frequency calculation with step size.").default_value(std::vector<double>{1, 1e-4}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-p", "--print").help("-- Printing options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-e", "--export").help("-- Export options.").default_value<std::vector<std::string>>({}).append();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{1, 1e-5}).nargs(2).scan<'g', double>();
    program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-o", "--optimize").help("-- Optimize the provided system.").default_value(1e-4).scan<'g', double>();

    // add positional arguments to the SCAN argument parser
    program.at<argparse::ArgumentParser>("scan").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").add_argument("-o", "--output").help("-- Output of the PES energies.").default_value("pes.dat");

    // add positional arguments to the SCAN HF argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();

    // add positional arguments to the SCAN MP2 argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the SCAN CI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the SCAN CIS argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the SCAN CID argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the SCAN CISD argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the SCAN FCI argument parser
    program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);

    // add positional arguments to the MD argument parser
    program.at<argparse::ArgumentParser>("md").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.5).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").add_argument("-o", "--output").help("-- Output of the trajectory.").default_value("trajectory.xyz");

    // add positional arguments to the MD HF argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-d", "--diis").help("-- Start iteration and history length for DIIS algorithm.").default_value(std::vector<int>{3, 5}).nargs(2).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-i", "--iters").help("-- Maximum number of iterations in SCF loop.").default_value(100).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-t", "--thresh").help("-- Threshold for conververgence in SCF loop.").default_value(1e-12).scan<'g', double>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD MP2 argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD CI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("excitations").help("-- Help message.").nargs(argparse::nargs_pattern::at_least_one).scan<'i', int>();
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD CIS argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD CID argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD CISD argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the MD FCI argument parser
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").add_argument("-g", "--gradient").help("-- Analytical (0) or numerical (1) gradient calculation with step size.").default_value(std::vector<double>{0, 1e-5}).nargs(2).scan<'g', double>();

    // add positional arguments to the QD argument parser
    program.at<argparse::ArgumentParser>("qd").add_argument("-h", "--help").help("-- Help message.").default_value(false).implicit_value(true);
    program.at<argparse::ArgumentParser>("qd").add_argument("-s", "--step").help("-- MD time step in atomic units.").default_value(0.1).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-i", "--iters").help("-- Number of iterations in dynamics.").default_value(1000).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-n", "--nstate").help("-- Number of states to consider.").default_value(3).scan<'i', int>();
    program.at<argparse::ArgumentParser>("qd").add_argument("-o", "--output").help("-- Output of the wavefunction.").default_value("wavefunction.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-f", "--potfile").help("-- File with the PES.").default_value("pes.dat");
    program.at<argparse::ArgumentParser>("qd").add_argument("-t", "--thresh").help("-- Threshold for conververgence in ITP loop.").default_value(1e-8).scan<'g', double>();
    program.at<argparse::ArgumentParser>("qd").add_argument("--no-real").help("-- Help message.").default_value(false).implicit_value(true);

    // parse the arguments
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
    for (const auto& parser : parsers) {
        try {if (parser.get<bool>("-h")) {
            std::cout << parser.help().str(); exit(EXIT_SUCCESS);
        }} catch (std::exception) {}
    }

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

    // catch unimplemented errors
    if (program.get<int>("-s") != 1) {
        if (program.is_subcommand_used("hf")) {
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) throw std::runtime_error("MP2 NOT IMPLEMENTED FOR UHF");
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) throw std::runtime_error("CI NOT IMPLEMENTED FOR UHF");
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) throw std::runtime_error("CIS NOT IMPLEMENTED FOR UHF");
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) throw std::runtime_error("CID NOT IMPLEMENTED FOR UHF");
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) throw std::runtime_error("CISD NOT IMPLEMENTED FOR UHF");
            if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) throw std::runtime_error("FCI NOT IMPLEMENTED FOR UHF");
        } else if (program.is_subcommand_used("scan")) {
            if (program.at<argparse::ArgumentParser>("scan").is_subcommand_used("hf")) {
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) throw std::runtime_error("MP2 NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) throw std::runtime_error("CI NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) throw std::runtime_error("CIS NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) throw std::runtime_error("CID NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) throw std::runtime_error("CISD NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) throw std::runtime_error("FCI NOT IMPLEMENTED FOR UHF");
            }
        } else if (program.is_subcommand_used("md")) {
            if (program.at<argparse::ArgumentParser>("md").is_subcommand_used("hf")) {
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) throw std::runtime_error("MP2 NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) throw std::runtime_error("CI NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) throw std::runtime_error("CIS NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) throw std::runtime_error("CID NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) throw std::runtime_error("CISD NOT IMPLEMENTED FOR UHF");
                if (program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) throw std::runtime_error("FCI NOT IMPLEMENTED FOR UHF");
            }
        }
    }
}

Distributor::~Distributor() {
    std::cout << "\n"; Printer::Title(std::string("TOTAL EXECUTION TIME: ") + Timer::Format(Timer::Elapsed(start)).c_str());
}

void Distributor::run() {
    // initialize the system
    std::string basis = program.get("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
    std::ifstream stream(program.get("-f")); system = System(stream, basis, program.get<int>("-c"), program.get<int>("-s")); stream.close();

    // print the initial stuff
    Printer::Initial(program, system);

    // calculate the integrals if needed
    if (program.is_subcommand_used("ints") || program.is_subcommand_used("hf")) integrals();

    // optimize the molecule
    if (program.is_subcommand_used("hf")) {
        if (program.at<argparse::ArgumentParser>("hf").is_used("-o")) {
            if (program.get<int>("-s") == 1) rhfo(program.at<argparse::ArgumentParser>("hf"));
            else uhfo(program.at<argparse::ArgumentParser>("hf"));
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").is_used("-o")) {
                if (program.get<int>("-s") == 1) rmp2o(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2"));
            }
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").is_used("-o")) {
                if (program.get<int>("-s") == 1) rcio(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci"));
            }
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").is_used("-o")) {
                if (program.get<int>("-s") == 1) rcio(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis"));
            }
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").is_used("-o")) {
                if (program.get<int>("-s") == 1) rcio(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid"));
            }
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").is_used("-o")) {
                if (program.get<int>("-s") == 1) rcio(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd"));
            }
        } else if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) {
            if (program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").is_used("-o")) {
                if (program.get<int>("-s") == 1) rcio(program.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci"));
            }
        }
    }

    // distribute the calculations
    if (program.is_subcommand_used("scan")) scan(program.at<argparse::ArgumentParser>("scan"));
    if (program.is_subcommand_used("md")) dynamics(program.at<argparse::ArgumentParser>("md"));
    if (program.is_subcommand_used("qd")) qdyn(program.at<argparse::ArgumentParser>("qd"));
    if (program.is_subcommand_used("hf")) {
        if (program.get<int>("-s") == 1) rhfrun(program.at<argparse::ArgumentParser>("hf"));
        else uhfrun(program.at<argparse::ArgumentParser>("hf"));
    }
}

void Distributor::rcirun(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the CI method header and define J in MO basis
    std::cout << "\n" + std::string(104, '-') + "\nRESTRICTED CONFIGURATION INTERACTION (";
    Matrix Hms; Tensor<4> Jms; CI::ResultsRestricted rcires;

    // print the method name
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) std::cout << "CISD)\n" << std::string(104, '-') + "\n";
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) std::cout << "CIS)\n" << std::string(104, '-') + "\n";
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) std::cout << "CID)\n" << std::string(104, '-') + "\n";
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) std::cout << "FCI)\n" << std::string(104, '-') + "\n";
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) std::cout << "CI)\n" << std::string(104, '-') + "\n";

    // transform the coulomb tensor
    std::cout << "\nCOULOMB INT IN MS BASIS: " << std::flush; TIME(Jms = Transform::CoulombSpin(system.ints.J, rhfres.C)) std::cout << " " << Eigen::MemTensor(Jms);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "jms")) {std::cout << "\n" << Jms << "\n";} std::cout << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "jms")) Eigen::Write("JMS.mat", Jms);
    
    // calculate the Hamiltonian in molecular sorbital basis
    std::cout << "HAMILTONIAN IN MS BASIS: " << std::flush; TIME(Hms = Transform::OneelecSpin(system.ints.T + system.ints.V, rhfres.C)) std::cout << " " << Eigen::MemMatrix(Hms);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "hms")) {std::cout << "\n" << Hms;} std::cout << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "hms")) Eigen::Write("HMS.mat", Hms);

    // do the calculation
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) rcires = CI({rhfres, parser.get<std::vector<int>>("excitations")}).rci(system, Hms, Jms);
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) rcires = CI({rhfres, {1, 2}}).rci(system, Hms, Jms);
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) rcires = CI({rhfres, {1}}).rci(system, Hms, Jms);
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) rcires = CI({rhfres, {2}}).rci(system, Hms, Jms);
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) rcires = CI({rhfres, {}}).rci(system, Hms, Jms);

    // print/save the result matrices
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "cih")) std::cout << "\nCI HAMILTONIAN\n" << rcires.H << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "cie")) std::cout << "\nCI ENERGIES\n" << Matrix(rcires.eig) << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "cic")) std::cout << "\nCI EXPANSION COEFFICIENTS\n" << rcires.C << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "cih")) Eigen::Write("CIH.mat", rcires.H);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "cie")) Eigen::Write("CIEPS.mat", Matrix(rcires.eig));
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "cic")) Eigen::Write("CIC.mat", rcires.C);

    // print the gradient and norm
    std::cout << "\nCI CORRELATION ENERGY: " << rcires.Ecorr << std::endl;
    std::cout << "FINAL CI ENERGY: " << rhfres.E + rcires.Ecorr << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) rcig(parser, rhfres);
    if (parser.is_used("-f")) rcif(parser, rhfres);
}

void Distributor::rcif(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for frequency calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED CI");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED CI");} Matrix H; CI::OptionsRestricted rciopt;

    // choose the correct method
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) rciopt.excits = parser.get<std::vector<int>>("excitations");
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) rciopt.excits = {1, 2};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) rciopt.excits = {1};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) rciopt.excits = {2};

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::ECI(rhfopt, rciopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RCI NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED CI FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rcig(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED CI");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED CI");} Matrix G; 

    // extract the RHF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb")); CI::OptionsRestricted rciopt;

    // choose the correct method
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) rciopt.excits = parser.get<std::vector<int>>("excitations");
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) rciopt.excits = {1, 2};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) rciopt.excits = {1};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) rciopt.excits = {2};

    // perform the CI gradient calculation
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::ECI(rhfopt, rciopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RCI NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rcio(argparse::ArgumentParser& parser) {
    // print the RHF optimization header
    std::cout << std::endl; Printer::Title("RESTRICTED CI OPTIMIZATION"); CI::OptionsRestricted rciopt;

    // extract the RHF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // choose the correct method
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) rciopt.excits = parser.get<std::vector<int>>("excitations");
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) rciopt.excits = {1, 2};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) rciopt.excits = {1};
    if (program.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) rciopt.excits = {2};

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGCI(rhfopt, rciopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::rhfrun(argparse::ArgumentParser& parser) {
    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK METHOD");
    std::printf("\n-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", rhfopt.maxiter, rhfopt.thresh, rhfopt.diis.start, rhfopt.diis.keep);

    // perform the RHF calculation
    auto rhfres = HF(rhfopt).rscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print/export the resulting RHF matrices and energies
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "eps")) std::cout << "\nORBITAL ENERGIES\n" << Matrix(rhfres.eps) << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "c")) std::cout << "\nCOEFFICIENT MATRIX\n" << rhfres.C << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "d")) std::cout << "\nDENSITY MATRIX\n" << rhfres.D << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "eps")) Eigen::Write("EPS.mat", Matrix(rhfres.eps));
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "c")) Eigen::Write("C.mat", rhfres.C);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "d")) Eigen::Write("D.mat", rhfres.D);
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << rhfres.E << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) rhfg(parser, rhfres);
    if (parser.is_used("-f")) rhff(parser, rhfres);

    // post RHF methods
    if (parser.is_subcommand_used("cisd")) rcirun(parser.at<argparse::ArgumentParser>("cisd"), rhfres);
    if (parser.is_subcommand_used("mp2")) rmp2run(parser.at<argparse::ArgumentParser>("mp2"), rhfres);
    if (parser.is_subcommand_used("cis")) rcirun(parser.at<argparse::ArgumentParser>("cis"), rhfres);
    if (parser.is_subcommand_used("cid")) rcirun(parser.at<argparse::ArgumentParser>("cid"), rhfres);
    if (parser.is_subcommand_used("fci")) rcirun(parser.at<argparse::ArgumentParser>("fci"), rhfres);
    if (parser.is_subcommand_used("ci")) rcirun(parser.at<argparse::ArgumentParser>("ci"), rhfres);
}

void Distributor::rhff(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED HARTREE-FOCK");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED HARTREE-FOCK");} Matrix H; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RHF NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rhfg(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED HARTREE-FOCK");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED HARTREE-FOCK");} Matrix G; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(rhfopt, rhfres.D));
    else G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, rhfres);

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rhfo(argparse::ArgumentParser& parser) {
    // print the RHF optimization header
    std::cout << std::endl; Printer::Title("RESTRICTED HARTREE-FOCK OPTIMIZATION");

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGHF(rhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::uhfrun(argparse::ArgumentParser& parser) {
    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // print the RHF method header
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK METHOD");
    std::printf("\n-- MAXITER: %d, THRESH: %.2e\n-- DIIS: [START: %d, KEEP: %d]\n", uhfopt.maxiter, uhfopt.thresh, uhfopt.diis.start, uhfopt.diis.keep);

    // perform the RHF calculation
    auto uhfres = HF(uhfopt).uscf(system, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));

    // print/save the resulting RHF matrices and energies
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "epsa")) std::cout << "\nORBITAL ENERGIES (ALPHA)\n" << Matrix(uhfres.epsa) << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "epsb")) std::cout << "\nORBITAL ENERGIES (BETA)\n" << Matrix(uhfres.epsb) << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "ca")) std::cout << "\nCOEFFICIENT MATRIX (ALPHA)\n" << uhfres.Ca << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "cb")) std::cout << "\nCOEFFICIENT MATRIX (BETA)\n" << uhfres.Cb << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "da")) std::cout << "\nDENSITY MATRIX (ALPHA)\n" << uhfres.Da << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "db")) std::cout << "\nDENSITY MATRIX (BETA)\n" << uhfres.Db << std::endl;
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "epsa")) Eigen::Write("EPSA.mat", Matrix(uhfres.epsa));
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "epsb")) Eigen::Write("EPSB.mat", Matrix(uhfres.epsb));
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "ca")) Eigen::Write("CA.mat", uhfres.Ca);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "cb")) Eigen::Write("CB.mat", uhfres.Cb);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "da")) Eigen::Write("DA.mat", uhfres.Da);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "db")) Eigen::Write("DB.mat", uhfres.Db);
    std::cout << "\nTOTAL NUCLEAR REPULSION ENERGY: " << Integral::Repulsion(system) << std::endl;
    std::cout << "FINAL HARTREE-FOCK ENERGY: " << uhfres.E << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) uhfg(parser, uhfres);
    if (parser.is_used("-f")) uhff(parser, uhfres);
}

void Distributor::uhff(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR UNRESTRICTED HARTREE-FOCK");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR UNRESTRICTED HARTREE-FOCK");} Matrix H; 

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR UHF NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::uhfg(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR UNRESTRICTED HARTREE-FOCK");} Matrix G; 

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // calculate the numerical or analytical gradient
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EHF(uhfopt, 0.5 * (uhfres.Da + uhfres.Db)));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR UHF NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::uhfo(argparse::ArgumentParser& parser) {
    // print the UHF optimization header
    std::cout << std::endl; Printer::Title("UNRESTRICTED HARTREE-FOCK OPTIMIZATION");

    // extract the HF options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGHF(uhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::rmp2run(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the MP2 method header and declare JMO tensor
    std::cout << std::endl; Printer::Title("RESTRICTED MØLLER-PLESSET PERTRUBATION THEORY"); Tensor<4> Jmo;

    // transform the coulomb tensor to MO basis
    std::cout << "\nCOULOMB TENSOR IN MO BASIS: " << std::flush; TIME(Jmo = Transform::Coulomb(system.ints.J, rhfres.C))
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "jmo")) {std::cout << "\n" << Jmo;} std::cout << "\n";
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-e"), "jmo")) Eigen::Write("JMO.mat", Jmo);

    // perform the MP2 calculation
    double Ecorr = MP({rhfres}).rmp2(system, Jmo);

    // print the gradient and it's norm
    std::cout << "\nMP2 CORRELATION ENERGY: " << Ecorr << std::endl << "FINAL ";
    std::cout << "MP2 ENERGY: " << rhfres.E + Ecorr << std::endl;

    // gradient and frequency
    if (parser.is_used("-g")) rmp2g(parser, rhfres);
    if (parser.is_used("-f")) rmp2f(parser, rhfres);
}

void Distributor::rmp2f(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for frequency calculation and define the hessian matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL FREQUENCIES FOR RESTRICTED MP2");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL FREQUENCIES FOR RESTRICTED MP2");} Matrix H; 

    // extract the HF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the Hessian calculation
    if (parser.get<std::vector<double>>("-f").at(0)) H = Hessian({parser.get<std::vector<double>>("-f").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL HESSIAN FOR RMP2 NOT IMPLEMENTED");

    // print the Hessian and it's norm
    Printer::Mat("\nNUCLEAR HESSIAN", H); std::printf("\nHESSIAN NORM: %.2e\n", H.norm());

    // calculate the frequencies
    Vector freq = Hessian({parser.get<std::vector<double>>("-f").at(1)}).frequency(system, H);

    // print the frequencies
    std::cout << std::endl; Printer::Title("RESTRICTED MP2 FREQUENCY ANALYSIS");
    Printer::Mat("\nVIBRATIONAL FREQUENCIES", freq);
}

void Distributor::rmp2g(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres) {
    // print the header for gradient calculation and define the gradient matrix
    if (parser.get<std::vector<double>>("-g").at(0)) {std::cout << std::endl; Printer::Title("NUMERICAL GRADIENT FOR RESTRICTED MP2");}
    else {std::cout << std::endl; Printer::Title("ANALYTICAL GRADIENT FOR RESTRICTED MP2");} Matrix G; 

    // extract the RHF options
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // perform the MP2 gradient calculation
    if (parser.get<std::vector<double>>("-g").at(0)) G = Gradient({parser.get<std::vector<double>>("-g").at(1)}).get(system, Lambda::EMP2(rhfopt, rhfres.D));
    else throw std::runtime_error("ANALYTICAL GRADIENT FOR RMP2 NOT IMPLEMENTED");

    // print the gradient and it's norm
    Printer::Mat("\nNUCLEAR GRADIENT", G); std::printf("\nGRADIENT NORM: %.2e\n", G.norm());
}

void Distributor::rmp2o(argparse::ArgumentParser& parser) {
    // extract the HF options and print the MP2 optimization method header
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    std::cout << std::endl; Printer::Title("RESTRICTED MP2 OPTIMIZATION");

    // perform the optimization
    system = Optimizer({parser.get<double>("-o")}).optimize(system, Lambda::EGMP2(rhfopt, parser.get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf())));

    // print the optimized coordinates and distances
    Printer::Mat("\nOPTIMIZED SYSTEM COORDINATES", system.coords);
    if (Utility::VectorContains<std::string>(parser.get<std::vector<std::string>>("-p"), "dist")) Printer::Mat("\nOPTIMIZED DISTANCE MATRIX", system.dists);
}

void Distributor::scan(argparse::ArgumentParser& parser) {
    // print the energy scan header
    std::cout << std::endl; Printer::Title("ENERGY SCAN"); std::printf("\nITER       Eel [Eh]           TIME    \n");

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // define the anonymous function for gradient
    std::function<double(System)> efunc; std::vector<double> energies;

    // get the energy and gradient function
    if (parser.is_subcommand_used("hf") && program.get<int>("-s") == 1) {
        efunc = Lambda::EHF(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            efunc = Lambda::EMP2(rhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) {
            efunc = Lambda::ECI(rhfopt, {{}, program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").get<std::vector<int>>("excitations")}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) {
            efunc = Lambda::ECI(rhfopt, {{}, {1}}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) {
            efunc = Lambda::ECI(rhfopt, {{}, {2}}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) {
            efunc = Lambda::ECI(rhfopt, {{}, {1, 2}}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) {
            efunc = Lambda::ECI(rhfopt, {{}, {}}, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.is_subcommand_used("hf") && program.get<int>("-s") != 1) {
        efunc = Lambda::EHF(uhfopt, Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
    }

    // count the number of geometries
    std::ifstream stream(program.get("-f")); std::string line; int lines;
    for (lines = 0; std::getline(stream, line); lines++);
    int geoms = lines / (system.atoms.size() + 2);
    stream.clear(), stream.seekg(0);

    // perform the scan
    for (int i = 0; i < geoms; i++) {
        // create the system
        std::string basis = program.get("-b"); std::replace(basis.begin(), basis.end(), '*', 's'), std::replace(basis.begin(), basis.end(), '+', 'p');
        system = System(stream, basis, program.get<int>("-c"), program.get<int>("-s")); Timer::Timepoint start = Timer::Now();

        // start the timer and calculate the energy
        double E = efunc(system); energies.push_back(E);

        // print the energy
        std::printf("%4d %20.14f %s\n", i + 1, E, Timer::Format(Timer::Elapsed(start)).c_str());
    }
    
    // save the file
    std::ofstream file(parser.get("-o"));
    file << std::fixed << std::setprecision(14) << "# i              E\n";
    for (size_t i = 0; i < energies.size(); i++) {
        file << std::setw(5) << i << " " << std::setw(20) << energies.at(i) << "\n";
    }
}

void Distributor::dynamics(argparse::ArgumentParser& parser) {
    // print the dynamics header
    std::cout << std::endl; Printer::Title("MOLECULAR DYNAMICS");

    // extract the HF method options
    auto uhfopt = HF::OptionsUnrestricted::Load(program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));
    auto rhfopt = HF::OptionsRestricted::Load(program.at<argparse::ArgumentParser>("md").at<argparse::ArgumentParser>("hf"), program.get<bool>("--no-coulomb"));

    // define the anonymous function for gradient
    std::function<std::tuple<double, Matrix>(System&)> egfunc;

    // get the energy and gradient function
    if (parser.is_subcommand_used("hf") && program.get<int>("-s") == 1) {
        egfunc = Lambda::EGHF(rhfopt, parser.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("mp2")) {
            egfunc = Lambda::EGMP2(rhfopt, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("mp2").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("ci")) {
            egfunc = Lambda::EGCI(rhfopt, {{}, program.at<argparse::ArgumentParser>("scan").at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").get<std::vector<int>>("excitations")}, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("ci").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cis")) {
            egfunc = Lambda::EGCI(rhfopt, {{}, {1}}, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cis").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cid")) {
            egfunc = Lambda::EGCI(rhfopt, {{}, {2}}, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cid").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("cisd")) {
            egfunc = Lambda::EGCI(rhfopt, {{}, {1, 2}}, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("cisd").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        } else if (parser.at<argparse::ArgumentParser>("hf").is_subcommand_used("fci")) {
            egfunc = Lambda::EGCI(rhfopt, {{}, {}}, parser.at<argparse::ArgumentParser>("hf").at<argparse::ArgumentParser>("fci").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
        }
    } else if (parser.is_subcommand_used("hf") && program.get<int>("-s") != 1) {
        egfunc = Lambda::EGHF(uhfopt, parser.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g"), Matrix::Zero(system.shells.nbf(), system.shells.nbf()));
    }

    // perform the dynamics
    MD({parser.get<int>("-i"), parser.get<double>("-s"), parser.get("-o")}).run(system, egfunc);
}

void Distributor::qdyn(argparse::ArgumentParser& parser) {
    // print the dynamics header
    std::cout << std::endl; Printer::Title("QUANTUM DYNAMICS");

    // perform the dynamics
    QD::Results qdres = QD({parser.get("-f"), parser.get<int>("-i"), parser.get<int>("-n"), parser.get<double>("-s"), parser.get<double>("-t"), parser.get<bool>("--no-real")}).run(system);

    // print the energies
    if (parser.get<bool>("--no-real")) std::cout << "\nIMAGINARY TIME PROPAGATION ENERGIES\n" << Matrix(qdres.energy) << std::endl;

    // save the wavefunction
    Utility::SaveWavefunction(parser.get<std::string>("-o"), qdres.r, qdres.states, qdres.energy);
}

void Distributor::integrals() {
    // print the integral calculation header
    std::cout << "\n"; Printer::Title("INTEGRAL CALCULATION");

    // calculate the overlap integral
    std::cout << "\nOVERLAP INTEGRAL: " << std::flush; TIME(system.ints.S = Integral::Overlap(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.S);
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "s")) std::cout << "\n" << system.ints.S << std::endl;
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-e"), "s")) Eigen::Write("S.mat", system.ints.S);
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-e"), "s")) Eigen::Write("S.mat", system.ints.S);

    // calculate the kinetic integral
    std::cout << "\nKINETIC INTEGRAL: " << std::flush; TIME(system.ints.T = Integral::Kinetic(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.T);
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "t")) std::cout << "\n" << system.ints.T << std::endl;
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-e"), "t")) Eigen::Write("T.mat", system.ints.T);
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-e"), "t")) Eigen::Write("T.mat", system.ints.T);

    // calculate the nuclear-electron attraction integral
    std::cout << "\nNUCLEAR INTEGRAL: " << std::flush; TIME(system.ints.V = Integral::Nuclear(system))
    std::cout << " " << Eigen::MemMatrix(system.ints.V);
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "v")) std::cout << "\n" << system.ints.V << std::endl;
    if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-e"), "v")) Eigen::Write("V.mat", system.ints.V);
    if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-e"), "v")) Eigen::Write("V.mat", system.ints.V);

    // calculate the electron-electron repulsion integral
    if (!program.get<bool>("--no-coulomb")) {std::cout << "\nCOULOMB INTEGRAL: " << std::flush; TIME(system.ints.J = Integral::Coulomb(system))}
    if (!program.get<bool>("--no-coulomb")) std::cout << " " << Eigen::MemTensor(system.ints.J);
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "j")) std::cout << "\n" << system.ints.J;
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-e"), "j")) Eigen::Write("J.mat", system.ints.J);
    if (!program.get<bool>("--no-coulomb") && program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-e"), "j")) Eigen::Write("J.mat", system.ints.J);

    // print new line
    std::cout << "\n";

    // if derivatives of the integrals are needed
    if (program.is_subcommand_used("ints") || ((program.at<argparse::ArgumentParser>("hf").is_used("-g") || program.at<argparse::ArgumentParser>("hf").is_used("-o")) && !program.at<argparse::ArgumentParser>("hf").get<std::vector<double>>("-g").at(0))) {
        // calculate the derivative of overlap integral
        std::cout << "\nFIRST DERIVATIVE OF OVERLAP INTEGRAL: " << std::flush; TIME(system.dints.dS = Integral::dOverlap(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dS);
        if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;
        if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "ds")) std::cout << "\n" << system.dints.dS << std::endl;

        // calculate the derivative of kinetic integral
        std::cout << "\nFIRST DERIVATIVE OF KINETIC INTEGRAL: " << std::flush; TIME(system.dints.dT = Integral::dKinetic(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dT);
        if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;
        if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dt")) std::cout << "\n" << system.dints.dT << std::endl;

        // calculate the derivative of nuclear-electron attraction integral
        std::cout << "\nFIRST DERIVATIVE OF NUCLEAR INTEGRAL: " << std::flush; TIME(system.dints.dV = Integral::dNuclear(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dV);
        if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;
        if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dv")) std::cout << "\n" << system.dints.dV << std::endl;

        // calculate the derivative of electron-electron repulsion integral
        std::cout << "\nFIRST DERIVATIVE OF COULOMB INTEGRAL: " << std::flush; TIME(system.dints.dJ = Integral::dCoulomb(system))
        std::cout << " " << Eigen::MemTensor(system.dints.dJ);
        if (program.is_subcommand_used("ints") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("ints").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;
        if (program.is_subcommand_used("hf") && Utility::VectorContains<std::string>(program.at<argparse::ArgumentParser>("hf").get<std::vector<std::string>>("-p"), "dj")) std::cout << "\n" << system.dints.dJ;

        // print new line
        std::cout << "\n";
    }
}
