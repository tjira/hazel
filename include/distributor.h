#pragma once

#include <sys/utsname.h>
#include <filesystem>

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
    // calculate integrals
    void integrals();

    // RHF distribution
    void rhfrun(argparse::ArgumentParser& parser); void rhff(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres);
    void rhfg(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres); void rhfo(argparse::ArgumentParser& parser);

    // RHF distribution
    void uhfrun(argparse::ArgumentParser& parser); void uhff(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres);
    void uhfg(argparse::ArgumentParser& parser, const HF::ResultsUnrestricted& uhfres); void uhfo(argparse::ArgumentParser& parser);

    // RCI distribution
    void rcirun(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres); void rcif(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres);
    void rcig(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres); void rcio(argparse::ArgumentParser& parser);

    // RMP2 distribution
    void rmp2run(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres); void rmp2f(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres);
    void rmp2g(argparse::ArgumentParser& parser, const HF::ResultsRestricted& rhfres); void rmp2o(argparse::ArgumentParser& parser);

    // MD and QD distribution
    void dynamics(argparse::ArgumentParser& parser); void qdyn(argparse::ArgumentParser& parser);

private:
    // argument parsers, timer and system
    std::vector<argparse::ArgumentParser> parsers;
    argparse::ArgumentParser program;
    Timer::Timepoint start;
    System system;
};
