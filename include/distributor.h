#pragma once

#include "integral.h"
#include "roothaan.h"
#include "timer.h"
#include "argparse.hpp"

class Distributor {

public:
    Distributor(int argc, char** argv);
    ~Distributor(); void run();

private:
    argparse::ArgumentParser program;
    Timer::Timepoint start;
};
