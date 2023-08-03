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

    // MD distribution
    void dynamics();

private:
    // printing options and parsers
    std::vector<std::string> print, ciprint, intsprint, hfprint, mdprint, mp2print;
    argparse::ArgumentParser program, ci, hf, ints, md, mdhf, mdmp2, mp2;

    // options and results of quantum methods
    HF::OptionsUnrestricted uhfopt; HF::ResultsUnrestricted uhfres;
    CI::OptionsRestricted rciopt; CI::ResultsRestricted rcires;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;

    // options and results of dynamics
    Dynamics::Options mdopt; Dynamics::Results mdres;

    // execution timestamp and system struct
    Timer::Timepoint start; System system;
};
