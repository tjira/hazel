#pragma once

#include <sys/utsname.h>
#include <filesystem>

#include "argparse.hpp"
#include "optimizer.h"
#include "dynamics.h"

#include <xc.h>

class Distributor {
public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    void integrals();

    // HF distribution
    void rhfrun(); void rhff();
    void rhfg(); void rhfo();

    // CI distribution
    void rcirun(); void rcif();
    void rcig(); void rcio();

    // MP2 distribution
    void rmp2run(); void rmp2f();
    void rmp2g(); void rmp2o();

    // MD distribution
    void dynamics();

private:
    // printing options and parsers
    argparse::ArgumentParser program, ints, hf, ci, md, mp2, mdhf, mdmp2;
    std::vector<std::string> print, ciprint, hfprint, mdprint, mp2print;

    // options and results of quantum methods
    CI::OptionsRestricted rciopt; CI::ResultsRestricted rcires;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;

    // options and results of dynamics
    Dynamics::Options mdopt; Dynamics::Results mdres;

    // execution timestamp and system struct
    Timer::Timepoint start; System system;
};
