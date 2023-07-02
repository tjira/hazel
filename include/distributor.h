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

private:
    argparse::ArgumentParser program, hf, ci, mp2;
    std::vector<std::string> print, save;
    Timer::Timepoint start;
};
