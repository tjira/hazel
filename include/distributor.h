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
    void hfrun(Data& data); void hff(Data& data);
    void hfg(Data& data); void hfo(Data& data);

    // MP2 distribution
    void mp2run(Data& data); void mp2f(Data& data);
    void mp2g(Data& data); void mp2o(Data& data);

private:
    std::vector<std::string> print, ciprint, hfprint, mp2print;
    argparse::ArgumentParser program, hf, ci, mp2;
    Timer::Timepoint start;

    Optimizer<HF>::Options opthfopt; Optimizer<HF>::Results opthfres;
    Optimizer<MP>::Options optmpopt; Optimizer<MP>::Results optmpres;
    Gradient<HF>::Options gradhfopt; Gradient<HF>::Results gradhfres;
    Gradient<MP>::Options gradmpopt; Gradient<MP>::Results gradmpres;
    Hessian<HF>::Options hesshfopt; Hessian<HF>::Results hesshfres;
    Hessian<MP>::Options hessmpopt; Hessian<MP>::Results hessmpres;
    CI::OptionsRestricted rciopt; CI::ResultsRestricted rcires;
    HF::OptionsRestricted rhfopt; HF::ResultsRestricted rhfres;
    MP::OptionsRestricted rmpopt; MP::ResultsRestricted rmpres;
    System system;
};
