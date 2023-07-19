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
    void hfrun(); void hff();
    void hfg(); void hfo();

    // MP2 distribution
    void mp2run(); void mp2f();
    void mp2g(); void mp2o();

private:
    std::vector<std::string> print, ciprint, hfprint, mp2print;
    argparse::ArgumentParser program, hf, ci, mp2;
    Timer::Timepoint start;

    Optimizer<HF>::OptionsRestricted opthfopt; Optimizer<MP>::OptionsRestricted optmpopt;
    Gradient<HF>::OptionsRestricted gradhfopt; Gradient<MP>::OptionsRestricted gradmpopt;
    Hessian<HF>::OptionsRestricted hesshfopt; Hessian<MP>::OptionsRestricted hessmpopt;
    CI::OptionsRestricted rciopt; CI::ResultsRestricted rcires;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;
    MP::OptionsRestricted rmpopt; MP::ResultsRestricted rmpres;
    System system;
};
