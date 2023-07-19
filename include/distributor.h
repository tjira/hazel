#pragma once

#include "argparse.hpp"
#include "optimizer.h"

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

private:
    // printing options and parsers
    std::vector<std::string> print, ciprint, hfprint, mp2print;
    argparse::ArgumentParser program, hf, ci, mp2;

    // options and results
    CI::OptionsRestricted rciopt; CI::ResultsRestricted rcires;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;

    // execution timestamp and system struct
    Timer::Timepoint start; System system;
};
