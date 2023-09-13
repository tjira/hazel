#pragma once

#include <sys/utsname.h>
#include <filesystem>

#include "argparse.hpp"
#include "optimizer.h"
#include "dynamics.h"
#include "lambda.h"
#include "qdyn.h"

#include <xc.h>

#ifndef CXXFLAGS
#define CXXFLAGS "---"
#endif

#ifndef DATADIR
#define DATADIR ""
#endif

class Distributor {
public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    // compute integrals
    void integrals();

    // RHF distribution
    void rhfrun(argparse::ArgumentParser& parser); void rhff(argparse::ArgumentParser& parser);
    void rhfg(argparse::ArgumentParser& parser); void rhfo(argparse::ArgumentParser& parser);

    // RHF distribution
    void uhfrun(argparse::ArgumentParser& parser); void uhff(argparse::ArgumentParser& parser);
    void uhfg(argparse::ArgumentParser& parser); void uhfo(argparse::ArgumentParser& parser);

    // RCI distribution
    void rcirun(argparse::ArgumentParser& parser); void rcif(argparse::ArgumentParser& parser);
    void rcig(argparse::ArgumentParser& parser); void rcio(argparse::ArgumentParser& parser);

    // RMP2 distribution
    void rmp2run(argparse::ArgumentParser& parser); void rmp2f(argparse::ArgumentParser& parser);
    void rmp2g(argparse::ArgumentParser& parser); void rmp2o(argparse::ArgumentParser& parser);

    // MD and QD distribution
    void dynamics(argparse::ArgumentParser& parser); void qdyn(argparse::ArgumentParser& parser);

private:
    // argument parsers, timer and system
    std::vector<argparse::ArgumentParser> parsers;
    argparse::ArgumentParser program;
    Timer::Timepoint start;
    System system;

    // options and results of quantum methods
    HF::OptionsUnrestricted uhfopt; HF::ResultsUnrestricted uhfres;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;
};
