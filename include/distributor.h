#pragma once

#include "integral.h"
#include "transform.h"
#include "mp.h"
#include "roothaan.h"
#include "timer.h"
#include "argparse.hpp"
#include <filesystem>

class Distributor {
public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    Integrals integrals(const System& system) const;

private:
    argparse::ArgumentParser program, hf, mp2;
    std::vector<std::string> print;
    Timer::Timepoint start;
};
