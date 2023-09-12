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
    void rhfrun(); void rhff();
    void rhfg(); void rhfo();

    // RHF distribution
    void uhfrun(); void uhff();
    void uhfg(); void uhfo();

    // RCI distribution
    void rcirun(); void rcif();
    void rcig(); void rcio();

    // RMP2 distribution
    void rmp2run(); void rmp2f();
    void rmp2g(); void rmp2o();

    // MD and QD distribution
    void dynamics(); void qdyn();

private:
    // printing options and parsers
    std::vector<std::string> print, ciprint, intsprint, hfprint, mdprint, mp2print, qdprint;
    argparse::ArgumentParser program, ints, hf, mp2, ci, qd, md, mdhf, mdmp2;

    // options and results of quantum methods
    HF::OptionsUnrestricted uhfopt; HF::ResultsUnrestricted uhfres;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;

    // execution timestamp and system struct
    Timer::Timepoint start; System system;
};
