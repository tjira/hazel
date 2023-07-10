#pragma once

#include "ci.h"
#include "mp.h"
#include "roothaan.h"
#include "transform.h"
#include "argparse.hpp"
#include <filesystem>

class Distributor {
public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    Data integrals(Data data) const;

    // Hartre-Fock distribution
    void hfrun(Data& data) const; void hff(Data& data) const;
    void hfg(Data& data) const; void hfo(Data& data) const;

    // MP2 distribution
    void mp2run(Data& data) const; void mp2f(Data& data) const;
    void mp2g(Data& data) const; void mp2o(Data& data) const;

private:
    std::vector<std::string> print, ciprint, hfprint, mp2print;
    argparse::ArgumentParser program, hf, ci, mp2;
    Timer::Timepoint start;
};
